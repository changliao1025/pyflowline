def check_if_duplicates(listOfElems):
    ''' Check if given list contains any duplicates '''    
    iFlag_unique = 1
    for elem in listOfElems:
        if listOfElems.count(elem) > 1:
            iFlag_unique = 0
            break
        else:
            pass
    
    return iFlag_unique