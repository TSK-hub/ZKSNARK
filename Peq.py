def Peq(Point1,Point2):
    if Point1['k']!=Point2['k']:
        raise ValueError('mismatched degree')
    if Point1['k']==Point2['k'] and Point1['x']==Point2['x'] and Point1['y']==Point2['y']:
        tf = 1
    else:
        tf = 0
    return tf
