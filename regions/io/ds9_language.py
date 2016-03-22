language = {'point': (coordinate, coordinate),
 'circle': (coordinate, coordinate, radius),
 'box': (coordinate, coordinate, width, height, angle),
}

circle(1.5, 3.6, 1.2)

def line_parser(line):
    region_type = 'circle' if 'circle' in line

    typer_parser(coordinate_string, language[region_type])

def type_parser(string, specfication):
    coord_list = []
    splitter = re.compile("[, ]")
    for element, element_parser in zip(splitter.split(string, specification):
        coord_list.append(element_parser(element))
