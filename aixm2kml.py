#!/usr/bin/env python3

import sys
import json
import math
import simplekml
import argparse
import warnings

from shapely.geometry import LineString, Point
from shapely.ops import split, nearest_points, snap
from bs4 import BeautifulSoup
from pyproj import Proj, transform
from pyproj import Transformer
import random
# Ignorer les avertissements RuntimeWarning de shapely
warnings.filterwarnings('ignore', category=RuntimeWarning, module='shapely')
warnings.filterwarnings("ignore", category=UserWarning, module='bs4')

def random_color():
    """Retourne une couleur aléatoire en format KML avec transparence."""
    a = 128  # Semi-transparent
    r = random.randint(0, 255)
    g = random.randint(0, 255)
    b = random.randint(0, 255)
    return simplekml.Color.rgb(a, r, g, b)

def geojson_to_kml(geojson, output_filename):
    kml = simplekml.Kml()
    for feature in geojson['features']:
        geom = feature['geometry']
        prop = feature['properties']
        name = prop.get('name', '')
        desc = prop.get('description', '')
        # Générez une couleur aléatoire pour ce volume
        color = random_color()
        if geom['type'] == 'Point':
            pnt = kml.newpoint(name=name, description=desc, coords=[tuple(geom['coordinates'])])
            pnt.style.iconstyle.color = color
        elif geom['type'] == 'LineString':
            lin = kml.newlinestring(name=name, description=desc, coords=geom['coordinates'])
            lin.style.linestyle.color = color
        elif geom['type'] == 'Polygon':
            pol = kml.newpolygon(name=name, description=desc, outerboundaryis=geom['coordinates'][0])
            pol.altitudemode = simplekml.AltitudeMode.absolute  # Définir le mode d'altitude à "absolute"
            pol.style.polystyle.color = color
        elif geom['type'] == 'Polyhedron':
            # Vous devrez peut-être ajuster cette partie pour gérer les polyèdres selon vos besoins
            pol = kml.newmultigeometry(name=name)
            for coords in geom['coordinates']:
                polygon = pol.newpolygon(outerboundaryis=coords)
                polygon.altitudemode = simplekml.AltitudeMode.absolute  # Définir le mode d'altitude à "absolute" pour chaque polygone dans le polyèdre
            pol.style.polystyle.color = color

    kml.save(output_filename)



def parse_args():
    parser = argparse.ArgumentParser(description="Convert AIXM data to KML.")
    parser.add_argument("--ahp", action="store_true", help="Export aerodromes/heliports")
    parser.add_argument("--obs", action="store_true", help="Export obstacles")
    parser.add_argument("--rcp", action="store_true", help="Export runway centers")
    parser.add_argument("--gbr", action="store_true", help="Export geographic borders")
    parser.add_argument("--abd", action="store_true", help="Export airspace boundaries")
    parser.add_argument("--uni", action="store_true", help="Export control towers")
    parser.add_argument("--gsd", action="store_true", help="Export gate stands")
    parser.add_argument("filename", help="Path to the AIXM file to process")
    return parser.parse_args()

def substring(geom, start_dist, end_dist, normalized=False):
    """Return a line segment between specified distances along a linear geometry.

    Negative distance values are taken as measured in the reverse
    direction from the end of the geometry. Out-of-range index
    values are handled by clamping them to the valid range of values.
    If the start distances equals the end distance, a point is being returned.
    If the normalized arg is True, the distance will be interpreted as a
    fraction of the geometry's length.

    from shapely 1.7
    """

    assert(isinstance(geom, LineString))
    
    # Filter out cases in which to return a point
    if start_dist == end_dist:
        return geom.interpolate(start_dist, normalized)
    elif not normalized and start_dist >= geom.length and end_dist >= geom.length:
        return geom.interpolate(geom.length, normalized)
    elif not normalized and -start_dist >= geom.length and -end_dist >= geom.length:
        return geom.interpolate(0, normalized)                    
    elif normalized and start_dist >= 1 and end_dist >= 1:
        return geom.interpolate(1, normalized)  
    elif normalized and -start_dist >= 1 and -end_dist >= 1:
        return geom.interpolate(0, normalized)

    start_point = geom.interpolate(start_dist, normalized)
    end_point = geom.interpolate(end_dist, normalized)
    
    min_dist = min(start_dist, end_dist)
    max_dist = max(start_dist, end_dist)
    if normalized:
        min_dist *= geom.length
        max_dist *= geom.length
    
    if start_dist < end_dist:
        vertex_list = [(start_point.x, start_point.y)]
    else:
        vertex_list = [(end_point.x, end_point.y)]
    coords = list(geom.coords)
    for i, p in enumerate(coords):
        pd = geom.project(Point(p))
        if min_dist < pd < max_dist:
            vertex_list.append(p)
        elif pd >= max_dist:
            break
    if start_dist < end_dist:
        vertex_list.append((end_point.x, end_point.y))
    else:
        vertex_list.append((start_point.x, start_point.y))
        # reverse direction result
        vertex_list = reversed(vertex_list)

    return LineString(vertex_list)


def geo2coordinates(o, latitude=None, longitude=None, recurse=True):
    if latitude:
        s = latitude
    else:
        s = o.find('geolat', recursive=recurse).string

    lat = s[:-1]
    if len(lat)==2 or lat[2]=='.': # DD[.dddd]
        lat = float(lat)
    elif len(lat)==4 or lat[4]=='.': # DDMM[.mmmm]
        lat = int(lat[0:2])+float(lat[2:])/60
    else: # DDMMSS[.sss]
        lat = int(lat[0:2])+int(lat[2:4])/60+float(lat[4:-1])/3600
    if s[-1] == 'S':
        lat = -lat

    if longitude:
        s = longitude
    else:
        s = o.find('geolong', recursive=recurse).string

    lon = s[:-1]
    if len(lon) == 3 or lon[3] == '.':
        lon = float(lon)
    elif len(lon) == 5 or lon[5] == '.':
        lon = int(lon[0:3])+float(lon[3:])/60
    else:
        lon = int(lon[0:3])+int(lon[3:5])/60+float(lon[5:-1])/3600
    if s[-1] == 'W':
        lon = -lon
    return([lon, lat])


def getfield(o, inputname, outputname=None):
    if outputname is None:
        outputname = inputname
    
    value = o.find(inputname.lower(), recursive=False)
    if value:
        return {outputname: value.string.replace('#','\n')}
    else:
        return None


def addfield(prop, field):
    if field:
        prop.update(field)
    return prop


def ahp2json(ahp):
    "Aerodrome / Heliport"

    # geometry
    geom = {"type": "Point", "coordinates": geo2coordinates(ahp)}

    # properties
    prop = dict()
    prop = addfield(prop, getfield(ahp, 'txtname', 'name'))
    prop = addfield(prop, getfield(ahp, 'codetype'))
    prop = addfield(prop, getfield(ahp, 'codeicao'))
    prop = addfield(prop, getfield(ahp, 'codeiata'))
    prop = addfield(prop, getfield(ahp, 'valelev','elevation'))
    prop = addfield(prop, getfield(ahp, 'uomdistver','vertical_unit'))
    prop = addfield(prop, getfield(ahp, 'txtdescrrefpt', 'description'))

    return {"type": "Feature", "geometry": geom, "properties": prop}


def obs2json(obs):
    "Obstacle"

    # geometry
    geom = {"type": "Point", "coordinates": geo2coordinates(obs)}

    # properties
    prop = dict()
    prop = addfield(prop, getfield(obs, 'txtname', 'name'))
    prop = addfield(prop, getfield(obs, 'txtdescrtype', 'description'))
    prop = addfield(prop, getfield(obs, 'txtDescrMarking', 'marked'))
    prop = addfield(prop, getfield(obs, 'codeLgt', 'light'))
    prop = addfield(prop, getfield(obs, 'valElev','elevation'))
    prop = addfield(prop, getfield(obs, 'valHgt','height'))
    prop = addfield(prop, getfield(obs, 'uomdistver','vertical_unit'))

    return {"type": "Feature", "geometry": geom, "properties": prop}

    # <Obs>
    #     <ObsUid mid="1577950">
    #         <geoLat>482936.00N</geoLat>
    #         <geoLong>0015526.00W</geoLong>
    #     </ObsUid>
    #     <txtName>22033</txtName>
    #     <txtDescrType>Pylône</txtDescrType>
    #     <codeGroup>N</codeGroup>
    #     <codeLgt>N</codeLgt>
    #     <txtDescrMarking>non balisé</txtDescrMarking>
    #     <codeDatum>U</codeDatum>
    #     <valElev>318</valElev>
    #     <valHgt>167</valHgt>
    #     <uomDistVer>FT</uomDistVer>
    # </Obs>


def rcp2json(o):
    "Runway Center line Position"

    # geometry
    geom = {"type": "Point", "coordinates": geo2coordinates(o.rcpuid)}

    # properties
    prop = dict()
    if o.ahpuid:
        prop = addfield(prop, getfield(o.ahpuid, 'codeid', 'codeicao'))
    if o.rwyuid:
        prop = addfield(prop, getfield(o.rwyuid, 'txtDesig', 'name'))
    prop = addfield(prop, getfield(o, 'valElev','elevation'))
    prop = addfield(prop, getfield(o, 'uomdistver','vertical_unit'))

    return {"type": "Feature", "geometry": geom, "properties": prop}

    # <Rcp>
    #     <RcpUid mid="1532169">
    #         <RwyUid mid="1528969">
    #             <AhpUid mid="1521170">
    #                 <codeId>LFMD</codeId>
    #             </AhpUid>
    #             <txtDesig>04/22</txtDesig>
    #         </RwyUid>
    #         <geoLat>433253.80N</geoLat>
    #         <geoLong>0065736.42E</geoLong>
    #     </RcpUid>
    #     <codeDatum>WGE</codeDatum>
    #     <valElev>10</valElev>
    #     <uomDistVer>FT</uomDistVer>
    # </Rcp>


def frange(start, stop, step):
    "Float range"
    if step > 0:
        while start < stop:
            yield start
            start += step
    else:
        while start > stop:
            yield start
            start += step


def xy2angle(x,y):
    if x == 0:
        if y>0:
            angle = pi/2
        else:
            angle = -pi/2
    else:
        angle = math.atan(y/x)
    if x<0:
        angle = angle + pi
    elif y<0:
        angle = (angle + 2 * pi)
    return angle % (2*pi)

def make_circle(lon, lat, radius, srs):
    g = []
    step = pi*2/120 if radius > 100 else pi*2/8
    
    transformer_to_srs = Transformer.from_proj(pWGS, srs)
    transformer_from_srs = Transformer.from_proj(srs, pWGS)
    
    center_x, center_y = transformer_to_srs.transform(lon, lat)
    for a in frange(0, pi*2, step):
        x = center_x + math.cos(a) * radius
        y = center_y + math.sin(a) * radius
        lon, lat = transformer_from_srs.transform(x, y)
        g.append([round(lon,6), round(lat,6)])
    return g


def abd2json(ase,gbr,o):
    "Airspace Border"
    # properties
    prop = dict()
    prop.update({"uid": o.abduid['mid']})
    if o.aseuid:
        prop = addfield(prop, getfield(o.aseuid, 'codetype'))
        prop = addfield(prop, getfield(o.aseuid, 'codeid'))

    if o.aseuid["mid"] in ase:
        a = ase[o.aseuid["mid"]]

        prop = addfield(prop, getfield(a, 'txtname', 'name'))
        prop = addfield(prop, getfield(a, 'codeclass', 'class'))
        prop = addfield(prop, getfield(a, 'codedistverupper', 'upper_type'))
        prop = addfield(prop, getfield(a, 'valdistverupper', 'upper_value'))
        prop = addfield(prop, getfield(a, 'uomdistverupper', 'upper_unit'))
        prop = addfield(prop, getfield(a, 'codedistverlower', 'lower_type'))
        prop = addfield(prop, getfield(a, 'valdistverlower', 'lower_value'))
        prop = addfield(prop, getfield(a, 'uomdistverlower', 'lower_unit'))
        prop = addfield(prop, getfield(a, 'txtrmk', 'remark'))
        # approximate altitudes in meters
        if a.uomdistverupper is not None:
            up = None
            if a.uomdistverupper.string == 'FL':
                up = float(a.valdistverupper.string) * ft * 100
            elif a.uomdistverupper.string == 'FT':
                up = float(a.valdistverupper.string) * ft
            elif a.uomdistverupper.string == 'M':
                up = float(a.valdistverupper.string)
            if up is not None:
                prop.update({"upper_m": int(up)})
        if a.uomdistverlower is not None:
            low = None
            if a.uomdistverlower.string == 'FL':
                low = float(a.valdistverlower.string) * ft * 100
            elif a.uomdistverlower.string == 'FT':
                low = float(a.valdistverlower.string) * ft
            elif a.uomdistverlower.string == 'M':
                low = float(a.valdistverlower.string)
            if low is not None:
                prop.update({"lower_m": int(low)})

    # geometry
    g = []
    if o.circle:
        lon_c, lat_c = geo2coordinates(o.circle, latitude=o.geolatcen.string, longitude=o.geolongcen.string)
        radius = float(o.valradius.string)
        if o.uomradius.string == 'NM':
            radius = radius * nm
        if o.uomradius.string == 'KM':
            radius = radius * 1000
        g = make_circle(lon_c, lat_c, radius, Proj(proj='ortho', lat_0=lat_c, lon_0=lon_c))

        # Si upper_value est différent de lower_value, alors c'est un volume
        if prop.get("upper_m") != prop.get("lower_m"):
            # Assumant que g contient maintenant le cercle 2D à l'altitude de base
            base_circle = [[lon, lat, prop.get("lower_m")] for lon, lat in g]
            top_circle = [[lon, lat, prop.get("upper_m")] for lon, lat in g]

            # Générer maintenant les côtés du cylindre en connectant les cercles de base et du haut
            sides = []
            for i in range(len(base_circle)):
                side = [
                    base_circle[i],
                    top_circle[i],
                    top_circle[(i + 1) % len(base_circle)],
                    base_circle[(i + 1) % len(base_circle)]
                ]
                sides.append(side)

            # Le polyèdre final a la base, le haut et les côtés
            polyhedron = [base_circle, top_circle] + sides

            geom = {"type": "Polyhedron", "coordinates": polyhedron}
        else:
            geom = {"type": "Polygon", "coordinates": [g]}

        prop = addfield(prop, getfield(o.valradius, 'radius'))
        prop = addfield(prop, getfield(o.uomradius, 'radius_unit'))
    else:
        avx_list = o.find_all('avx')
        for avx_cur in range(0,len(avx_list)):
            avx = avx_list[avx_cur]
            codetype = avx.codetype.string
            if codetype in ['GRC', 'RHL']:
                # great-circle segment
                g.append(geo2coordinates(avx))
            elif codetype in ['CCA', 'CWA']:
                # arcs
                start = geo2coordinates(avx, recurse=False)
                if avx_cur+1 == len(avx_list):
                    stop = g[0]
                else:
                    stop = geo2coordinates(avx_list[avx_cur+1], recurse=False)
                center = geo2coordinates(avx,
                                         latitude=avx.geolatarc.string,
                                         longitude=avx.geolongarc.string)
                g.append(start)
                # convert to local meters
                srs = Proj(proj='ortho', lat_0=center[1], lon_0=center[0])
                # Création du transformer
                transformer_to_srs = Transformer.from_proj(pWGS, srs)
                # Transformation des coordonnées
                start_x, start_y = transformer_to_srs.transform(start[0], start[1])
                stop_x, stop_y = transformer_to_srs.transform(stop[0], stop[1])
                center_x, center_y = transformer_to_srs.transform(center[0], center[1])
                # recompute radius from center/start coordinates in local projection
                radius = math.sqrt(start_x**2+start_y**2)
                # start / stop angles
                start_angle = round(xy2angle(start_x-center_x, start_y-center_y),6)
                stop_angle = round(xy2angle(stop_x-center_x, stop_y-center_y),6)
                step = -0.025 if codetype == 'CWA' else 0.025
                if codetype == 'CWA' and stop_angle > start_angle:
                    stop_angle = stop_angle - 2*pi
                if codetype == 'CCA' and stop_angle < start_angle:
                    start_angle = start_angle - 2*pi
                for a in frange(start_angle+step/2, stop_angle-step/2, step):
                    x = center_x + math.cos(a) * radius
                    y = center_y + math.sin(a) * radius
                    # Création du transformer pour la transformation inverse
                    transformer_from_srs = Transformer.from_proj(srs, pWGS)
                    # Transformation des coordonnées
                    lon, lat = transformer_from_srs.transform(x, y)
                    g.append([lon, lat])
            elif codetype == 'FNT':
                # geographic borders
                start = geo2coordinates(avx)
                if avx_cur+1 == len(avx_list):
                    stop = g[0]
                else:
                    stop = geo2coordinates(avx_list[avx_cur+1])
                if avx.gbruid["mid"] in gbr:
                    fnt = gbr[avx.gbruid["mid"]]
                    start_d = fnt.project(Point(start[0], start[1]), normalized=True)
                    stop_d = fnt.project(Point(stop[0], stop[1]), normalized=True)
                    geom = substring(fnt, start_d, stop_d, normalized=True)
                    for c in geom.coords:
                        lon, lat = c
                        g.append([lon, lat])
                else:
                    print('!!! missing GBR', avx.gbruid["mid"])
                    g.append(start)
            else:
                g.append(geo2coordinates(avx))

        if (len(g)==0):
            print(o.prettify())
            geom = None
        elif len(g) == 1:
            geom = {"type": "Point", "coordinates": g[0]}
        elif len(g) == 2:
            geom = {"type": "LineString", "coordinates": g}
        else:
            base_polygon = [[lon, lat, prop.get("lower_m")] for lon, lat in g]
            top_polygon = [[lon, lat, prop.get("upper_m")] for lon, lat in g]
            # Si upper_value est différent de lower_value, alors c'est un volume
            if prop.get("upper_m") != prop.get("lower_m"):
                sides = []
                for i in range(len(base_polygon)):
                    side = [
                        base_polygon[i],
                        top_polygon[i],
                        top_polygon[(i + 1) % len(base_polygon)],
                        base_polygon[(i + 1) % len(base_polygon)]
                    ]
                    sides.append(side)
                # Le polyèdre final a la base, le sommet et les côtés
                polyhedron = [base_polygon, top_polygon] + sides             
                geom = {"type": "Polyhedron", "coordinates": polyhedron}
              
            else:
                base_polygon.append(base_polygon[0])  # Fermez le polygone
                geom = {"type": "Polygon", "coordinates": [base_polygon]}
      
    return {"type": "Feature", "geometry": geom, "properties": prop}

    # <Abd>
    #     <AbdUid mid="1570252">
    #         <AseUid mid="1562867">
    #             <codeType>SECTOR</codeType>
    #             <codeId>LFRRVU</codeId>
    #         </AseUid>
    #     </AbdUid>
    #     <Avx>
    #         <codeType>GRC</codeType>
    #         <geoLat>493400.00N</geoLat>
    #         <geoLong>0042700.00W</geoLong>
    #         <codeDatum>WGE</codeDatum>
    #     </Avx>
    #     <Avx>
    #         <codeType>GRC</codeType>
    #         <geoLat>500000.00N</geoLat>
    #         <geoLong>0020000.00W</geoLong>
    #         <codeDatum>WGE</codeDatum>
    #     </Avx>
    #     <Avx>
    #         <codeType>GRC</codeType>
    #         <geoLat>494343.00N</geoLat>
    #         <geoLong>0015825.00W</geoLong>
    #         <codeDatum>WGE</codeDatum>
    #     </Avx>
    #     <Avx>
    #         <codeType>GRC</codeType>
    #         <geoLat>485404.00N</geoLat>
    #         <geoLong>0024822.00W</geoLong>
    #         <codeDatum>WGE</codeDatum>
    #     </Avx>
    #     <Avx>
    #         <codeType>GRC</codeType>
    #         <geoLat>484607.00N</geoLat>
    #         <geoLong>0025950.00W</geoLong>
    #         <codeDatum>WGE</codeDatum>
    #     </Avx>
    #     <Avx>
    #         <codeType>GRC</codeType>
    #         <geoLat>484746.00N</geoLat>
    #         <geoLong>0040157.00W</geoLong>
    #         <codeDatum>WGE</codeDatum>
    #     </Avx>
    # </Abd>


def gbr2json(o):
    "Geographic borders"
    # geometry
    g = []
    l = []
    for gbv in o.find_all('gbv'):
        if gbv.codetype.string not in ['GRC', 'END']:
            print(gbv)
        g.append(geo2coordinates(gbv))
        l.append((g[-1][0], g[-1][1]))
    geom = { "type":"LineString", "coordinates": g }

    # properties
    prop = dict()
    prop = addfield(prop, getfield(o.gbruid, 'txtname', 'name'))
    prop = addfield(prop, getfield(o, 'codetype', 'type'))
    prop = addfield(prop, getfield(o, 'txtrmk', 'remark'))
    
    return ({"type": "Feature", "geometry": geom, "properties": prop} , l)

    # <Gbr>
    #     <GbrUid mid="1544998">
    #         <txtName>ZR:LE</txtName>
    #     </GbrUid>
    #     <codeType>OTHER</codeType>
    #     <Gbv>
    #         <codeType>GRC</codeType>
    #         <geoLat>424729.76N</geoLat>
    #         <geoLong>0000031.32W</geoLong>
    #         <codeDatum>WGE</codeDatum>
    #     </Gbv>
    #     ...


def tower2json(o):
    "Control towers"
    if o.codetype.string == 'TWR' and o.find('geolat'):
        # geometry
        geom = { "type":"Point", "coordinates": geo2coordinates(o) }

        # properties
        prop = dict()
        prop = addfield(prop, getfield(o.uniuid, 'txtname', 'name'))

        return {"type": "Feature", "geometry": geom, "properties": prop}
    else:
        print("!!! missing TWR coordinates", o)

    # <Uni>
    #     <UniUid mid="1524684">
    #         <txtName>LFBR MURET</txtName>
    #     </UniUid>
    #     <OrgUid mid="1520800">
    #         <txtName>FRANCE</txtName>
    #     </OrgUid>
    #     <AhpUid mid="1521106">
    #         <codeId>LFBR</codeId>
    #     </AhpUid>
    #     <codeType>TWR</codeType>
    #     <codeClass>OTHER</codeClass>
    #     <geoLat>432656.96N</geoLat>
    #     <geoLong>0011549.30E</geoLong>
    #     <codeDatum>WGE</codeDatum>
    # </Uni>
 
def gsd2json(o):
    "Gate stands"
    # geometry
    geom = { "type":"Point", "coordinates": geo2coordinates(o) }

    # properties
    prop = dict()
    prop = addfield(prop, getfield(o.gsduid, 'txtdesig', 'ref'))
    prop = addfield(prop, getfield(o.gsduid.apnuid.ahpuid, 'codeid', 'airport_ref'))

    return {"type": "Feature", "geometry": geom, "properties": prop}

    # <Gsd>
    #     <GsdUid mid="1583952">
    #         <ApnUid mid="1574122">
    #             <AhpUid mid="1520960">
    #                 <codeId>LFMN</codeId>
    #             </AhpUid>
    #             <txtName>LFMN-APRON</txtName>
    #         </ApnUid>
    #         <txtDesig>Y2</txtDesig>
    #     </GsdUid>
    #     <codeType>OTHER</codeType>
    #     <geoLat>433923.85N</geoLat>
    #     <geoLong>0071214.03E</geoLong>
    #     <codeDatum>WGE</codeDatum>
    # </Gsd>

pLocal = Proj('epsg:2154')
pWGS = Proj('epsg:4326')

#pLocal = Proj(init='epsg:2154')
#pWGS = Proj(init='epsg:4326')

nm = 1852   # Nautic Mile to meters
ft = 0.3048 # foot in meter
pi = 3.1415926


def extract_ahp(aixm):
    #print("extract ahp - aerodromes/heliports")
    out = []
    for o in aixm.find_all('ahp'):
        out.append(ahp2json(o))
    geojson_data = {"type":"FeatureCollection", "features": out}
    
    geojson_to_kml(geojson_data, 'aerodromes.kml')

def extract_obs(aixm):
    #print("extract obs - obstacles")
    out = []
    for o in aixm.find_all('obs'):
        out.append(obs2json(o))
    geojson_data = {"type":"FeatureCollection", "features": out}
    geojson_to_kml(geojson_data, 'obstacles.kml')

def extract_rcp(aixm):
    #print("extract rcp - runway centers")
    out = []
    for o in aixm.find_all('rcp'):
        out.append(rcp2json(o))
    geojson_data = {"type":"FeatureCollection", "features": out}
    geojson_to_kml(geojson_data, 'runway_center.kml')

def extract_ase(aixm):
    print("extract ase - airspace")
    ase = dict()
    for o in aixm.find_all('ase'):
        ase[o.aseuid['mid']] = o
    return ase

def extract_gbr(aixm):
    #print("extract gbr - geographic borders")
    gbr = dict()
    out = []
    for o in aixm.find_all('gbr'):
        j,l = gbr2json(o)
        out.append(j)
        gbr[o.gbruid['mid']] = LineString(l)
    geojson_data = {"type":"FeatureCollection", "features": out}
    geojson_to_kml(geojson_data, 'border.kml')
    return gbr

import os

def extract_abd(aixm):
    print("extract abd - airspace boundaries")
    ase = extract_ase(aixm)
    gbr = extract_gbr(aixm)

    # Assurez-vous que le répertoire "airspace" existe
    if not os.path.exists('airspace'):
        os.makedirs('airspace')

    for o in aixm.find_all('abd'):
        boundary = abd2json(ase, gbr, o)
        
        # Utilisez le code ou le nom pour nommer le fichier
        filename = boundary['properties'].get('codeid', 'unknown')
        
        # Évitez les caractères spéciaux dans le nom de fichier
        filename = ''.join(e for e in filename if e.isalnum())

        # Obtenez le type de géométrie
        geom_type = boundary['geometry']['type']

        # Assurez-vous que le sous-répertoire pour ce type de géométrie existe
        geom_dir = os.path.join('airspace', geom_type)
        if not os.path.exists(geom_dir):
            os.makedirs(geom_dir)
        
        # Créez le chemin complet du fichier avec le sous-répertoire
        kml_file_path = os.path.join(geom_dir, f"{filename}.kml")
        
        # Convertissez la limite en GeoJSON et sauvegardez-la en KML
        geojson_data = {"type": "FeatureCollection", "features": [boundary]}
        # Si le codeid est "LFTR200E", affichez geojson_data
        if filename == "LFTR200E":
            print("**************************************")
            print(geojson_data)   
        geojson_to_kml(geojson_data, kml_file_path)



def extract_uni(aixm):
    #print("extract uni - control towers")
    out = []
    for o in aixm.find_all('uni'):
        twr = tower2json(o)
        if twr:
            out.append(twr)
    geojson_data = {"type":"FeatureCollection", "features": out}
    geojson_to_kml(geojson_data, 'tower.kml')

def extract_gsd(aixm):
    #print("extract gsd - gate stands")
    out = []
    for o in aixm.find_all('gsd'):
        out.append(gsd2json(o))
    geojson_data = {"type":"FeatureCollection", "features": out}
    geojson_to_kml(geojson_data, 'gate_stand.kml')

if __name__ == "__main__":
    args = parse_args()
  
    
    print(f"Parsing XML file: {args.filename}")
    
    with open(args.filename, 'r') as f:
        aixm = BeautifulSoup(f, 'lxml')
    
    if args.ahp:
        print("extract ahp - aerodromes/heliports")
        extract_ahp(aixm)

    if args.obs:
        print("extract obs - obstacles")
        extract_obs(aixm)

    if args.rcp:
        print("extract rcp - runway centers")
        extract_rcp(aixm)

    if args.gbr:
        print("extract gbr - geographic borders")
        extract_gbr(aixm)

    if args.abd:
        print("extract abd - airspace boundaries")
        extract_abd(aixm)

    if args.uni:
        print("extract uni - control towers")
        extract_uni(aixm)

    if args.gsd:
        print("extract gsd - gate stands")
        extract_gsd(aixm)

    print("done")



