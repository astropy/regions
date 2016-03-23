import re
import string
import cgi
from ..shapes import circle, ellipse, rectangle
from ..core.pixcoord import PixCoord
import astropy.coordinates as coords



def parse_ds9(filename=None,save_comments=False):
    """
    parse a ds9 regions file and return a regions object

        example somple file contents:

        # Region file format: DS9 version 4.1
        # Filename: /Users/sosey/test_images/iacs01t4q_flt.fits[SCI]
        global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
        physical
        circle(622,605,20)
        circle(656,534,20)
        circle(528,555,20)
        circle(484,606,20)
        ellipse(267,601,55,24,0)
        ellipse(409,498,40,20,0) # color=magenta


    """
    if not filename:
        print("Please provide filename for ds9 region file")
        return ValueError

    shapes=["circle","ellipse","point","box"]
    special=["global"]
    wcs=["physical","fk5"]
    match_params = re.compile("([a-zA-Z]+)(=)([^=]+[$| }*\"*])+")

    newline=";"
    comments=list()

    with open(filename,'r') as fh:
        lines=fh.read().splitlines()

    #list of region objects which are returned
    all_regions=list()
    region_shape=None
    sprinkles=None

    #doesn't deal with quotes in parameter specs yet
    #this looping needs some serious work, I think a
    #combination of more regex and string expressions
    #would be a decent compromise for speed and readability
    for line in lines:
        if line.startswith("#"): #comment, save for later
            comments.append(line)
        else:
            this_region=None
            vertex=None
            param_list=list()
            if not line.startswith(string.punctuation):
                if newline in line:
                    morelines=line.split(newline)
                    #sketchy growth of list
                    for n in morelines:
                        lines.append(n.strip())
                else:
                    if line.count("#") == 1:
                        region_shape, sprinkles =line.split("#")
                        param_list=match_params.findall(sprinkles)
                    else:
                        region_shape=line

                    region=None
                    for shape in shapes:
                        if shape in region_shape:
                            region=shape
                            loc=region_shape.find("(")
                            if loc:
                                loc2=region_shape.find(")")
                                if "," in region_shape:
                                    splitter=","
                                else:
                                    splitter=" "
                                contents=region_shape[loc+1:loc2].split(splitter)
                                if ":" in region_shape:
                                    vertex=[coords.SkyCoord(contents[0],contents[1],unit="deg")]
                                else:
                                    vertex=[float(num) for num in contents]
                    if vertex:
                        if region is "circle":
                            if len(vertex) < 2:
                                vertex.append(contents[-1])
                                this_region=circle.CirclePixelRegion(vertex[0],vertex[1],params=param_list)
                            else:
                                x,y,radius=zip(vertex)
                                this_region=circle.CirclePixelRegion((x,y),radius,params=param_list)
                        if region is "ellipse":
                            this_region=ellipse.EllipsePixelRegion(vertex,params=param_list)
                        if region is "point":
                            x,y=zip(vertex)
                            this_region=PixCoord(x,y,params=param_list)
                        if region is "box":
                            this_region=rectangle.RectanglePixelRegion(vertex,params=param_list)

                if this_region:
                    all_regions.append(this_region) #(region,vertex,param_list)

    if (save_comments):
        while comments:
            all_regions.append(comments.pop())

    return all_regions
