import os

def get_package_data():
    data = [os.path.join('data', '*.fits'),
            os.path.join('baseline', '*.fits')]
    return {'regions.shapes.tests': data}
