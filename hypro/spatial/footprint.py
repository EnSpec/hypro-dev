from shapely.geometry import LineString, Polygon, MultiPolygon
from shapely.geometry import JOIN_STYLE

def footprint(igm, ps):
    
    # Construct linestring from edge points of IGM
    ordered_points = []
    ordered_points += [*igm[:2,:,0].T] # Port edge?
    ordered_points += [*igm[:2,-1,1:-1].T] # Ending edge
    ordered_points += [*igm[:2,::-1,-1].T] # Starbord edge?
    ordered_points += [*igm[:2,0,-2::-1].T] # Starting edge
    
    ### TODO: Filter out any NaN values!
    
    # Buffer the linestring by half the pixel size
    buffered_outline = LineString(ordered_points).buffer(0.5*abs(ps), join_style=JOIN_STYLE.bevel)
    # buffered_outline = LineString(ordered_points).buffer(0.5*ps, join_style=JOIN_STYLE.mitre, mitre_limit=ps)
    
    if type(buffered_outline) is Polygon:
        # Isolate the exterior geometry
        footprint = Polygon(buffered_outline.exterior)
    elif type(buffered_outline) is MultiPolygon:
        raise NotImplementedError
    
    return footprint