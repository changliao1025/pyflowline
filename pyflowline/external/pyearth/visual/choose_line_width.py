def choose_line_width(iFigure_width, iDPI):
    # Calculate line width based on the figure width and iDPI
    dummy = 0.01 * iFigure_width * iDPI / 80
    dummy1 = max(2.5, dummy)
    iWidth_out = min(1.5, dummy1)
    return iWidth_out

if __name__ == '__main__':
    iFigure_width = 8 #inch
    iDPI = 300 #dot per inch
    iWidth_out = choose_line_width(iFigure_width, iDPI)
    print(iWidth_out)