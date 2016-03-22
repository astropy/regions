import re
from regions.shapes.circle import CirclePixelRegion
from regions.shapes.ellipse import EllipsePixelRegion


def parse_ds9(filename=None):
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



    More specifically, region specifications consist of one or more lines containing:
      # comment until end of line
      global   keyword=value keyword=value  ... # set global value(s)
      # include the following file in the region descriptor
      @file
      # use the FITS image as a mask (cannot be used with other regions)
      @fitsimage
      # each region expression contains shapes separated by operators
      [region_expression1], [region_expression2], ...
      [region_expression], [region_expression], ...
    A single region expression consists of:

      # parens and commas are optional, as is the + sign
      [+-]shape(num , num , ...) OP1 shape num num num OP2 shape ...

    e.g.:

      ([+-]shape(num , num , ...) && shape num  num || shape(num, num)
      # a comment can come after a region -- reserved for local properties
      [+-]shape(num , num , ...)  # local properties go here, e.g. color=red
    Thus, a region descriptor consists of one or more region expressions or regions, separated by comas, new-lines, or semi-colons. Each region consists of one or more geometric shapes combined using standard boolean operation. Several types of shapes are supported, including:


    Syntax
    ------
    Region arguments may be separated with either a comma or space. Optional parentheses may be used a the beginning and end of a description.

    circle 100 100 10
    circle(100 100 10)
    circle(100,100,10)

    Comments
    --------
    All lines that begin with # are comments and will be ignored.

    # This is a comment

    Delimiter
    ---------
    All lines may be delimited with either a new-line or semi-colon.

    circle 100 100 10
    ellipse 200 200 20 40 ; box 300 300 20 40

    Header
    ------

    A DS9 region file may start with the following optional header:

    # Region file format: DS9 version 4.0

    Global Properties
    -----------------
    Global properties affect all regions unless a local property is specified. The global keyword is first, followed by a list of keyword = value pairs. Multiple global property lines may be used within a region file.

    global color=green font="helvetica 10 normal roman" edit=1 move=1 delete=1 highlite=1 include=1 wcs=wcs

    Local Properties
    ----------------
    Local properties start with a # after a region description and only affect the region it is specified with.

    physical;circle(504,513,20) # color=red text={This is a Circle}


    The arguments to region shapes can be floats or integers describing positions and sizes. They can be specified as pure numbers or using explicit formatting directives:

    position arguments

    [num]                   # context-dependent (see below)
    [num]d                  # degrees
    [num]r                  # radians
    [num]p                  # physical pixels
    [num]i                  # image pixels
    [num]:[num]:[num]       # hms for 'odd' position arguments
    [num]:[num]:[num]       # dms for 'even' position arguments
    [num]h[num]m[num]s      # explicit hms
    [num]d[num]m[num]s      # explicit dms
    size arguments

    [num]                   # context-dependent (see below)
    [num]"                  # arc sec
    [num]'                  # arc min
    [num]d                  # degrees
    [num]r                  # radians
    [num]p                  # physical pixels
    [num]i                  # image pixels

    """
    if not filename:
        print("Please provide filename for ds9 region file")
        return ValueError

    shapes=["circle","ellipse"]
    special=["global"]

    #seperation tokens
    tokens=["#","(",")","="]
    newline=";"
    comments=[]
    sline=[]

    #match hex [0-9A-Fa-f]
    #match letters [a-zA-Z]

    with open(filename,'r') as fh:
        lines=fh.read().splitlines()

    #list of region objects which are returned
    all_regions=list()
    this_region=None

    #doesn't deal with quotes in parameter specs yet
    for line in lines:
        this_region=None
        vertex=None
        param_list=list()
        span=re.search("[a-zA-Z]",line)
        if span.span()[0] == 0:#character found in first position
            if newline in line:
                sline=line.split(newline)
                #sketchy growth of list
                for newline in line:
                        lines.append(newline)
            else:
                sline=line.split()

                #find parameters
                for spec in sline:
                    if "=" in spec:
                        param,val=spec.split("=")
                        param_list.append((param,val))
                    else:
                        for key in special:
                            if key in spec:
                                region=key
                        else:
                            for shape in shapes:
                                if shape in spec:
                                    region=shape
                                    loc=spec.find("(")
                                    if loc:
                                        loc2=spec.find(")")
                                        vertex=[float(num) for num in spec[loc+1:loc2].split(",")]
                                        if region == "circle":
                                            x,y,radius=zip(vertex)
                                            this_region=CirclePixelRegion((x,y),radius)
                                        if region == "ellipse":
                                            this_region=EllipsePixelRegion(vertex)

                if this_region:
                    all_regions.append(this_region) #(region,vertex,param_list)

        elif tokens[0] in line[0]: #comment, save for later
            comments.append(line)

    while comments:
            all_regions.append(comments.pop())

    return all_regions
