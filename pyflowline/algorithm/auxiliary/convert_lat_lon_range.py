def convert_360_to_180(dLongitude_in):
    '''
    This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html
    '''
    a = int(dLongitude_in /180)
    dLongitude_out = dLongitude_in - a*360.0

    return dLongitude_out



def convert_180_to_360(dLongitude_in):
    '''
    This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html
    '''
    dLongitude_out = (dLongitude_in + 360.0) % 360.0

    return dLongitude_out