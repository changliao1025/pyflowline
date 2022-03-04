import os
import numpy as np
def text_reader_string( sFilename_in,\
     ncolumn_in = None, \
     nrow_in = None, \
     cDelimiter_in = None, \
     iFlag_remove_quota = None, \
     iSkipline_in =  None ):
    """
    Read a text based file
    sFilename_in,
    ncolumn_in = None, 
    nrow_in = None, 
    cDelimiter_in = None, 
    iSkipline_in =  None
    """

    aData_out = -1
    if os.path.isfile(sFilename_in):

        if ncolumn_in is not None:
            iFlag_column = 1
            ncolumn_out = ncolumn_in
        else:
            iFlag_column = 0
            
        if nrow_in is not None :
            #there is nothing
            nrow_out = nrow_in
            #print(' ')
        else :
            #get the total line number
            ifs = open(sFilename_in, "r")
            nrow_out = len(ifs.readlines())
            ifs.close()

        sLine=' '
        ifs = open(sFilename_in, "r")

        if iSkipline_in is not None:
            nrow_out =  nrow_out - iSkipline_in
            for i in range(iSkipline_in):
                ifs.readline()
        else:
            pass

        if iFlag_remove_quota is not None:
            iFlag_iFlag_remove_quota = 1
        else:
            iFlag_iFlag_remove_quota = 0

        #get delimiter
        if cDelimiter_in is not None:
            iFlag_delimiter = 1
        else:
            iFlag_delimiter = 0
            #cDelimiter_in = ' '
            pass
        
        
        if iFlag_delimiter == 1:
            if iFlag_column == 1:
                if ncolumn_out < 1 :
                    print( 'The file has no data!' )
                    pass

                aData_out = np.full( (nrow_out, ncolumn_out), '  ' , dtype=object )    
                for iRow in range(nrow_out):
                    sLine=(ifs.readline()).rstrip()
                    if iFlag_iFlag_remove_quota ==1:
                        sLine.replace('"','')
                    else:
                        pass    

                    aDummy = sLine.split(cDelimiter_in)
                    iDummy = len(aDummy)
                    if iDummy == ncolumn_out:
                        aData_out[iRow] = aDummy
                    else:
                        #something is missing
                        aDummy2  =np.full(ncolumn_out, '  ' , dtype=object)
                        aDummy2[0: iDummy]= aDummy
                        aData_out[iRow]= aDummy2
                        pass

            else :
                sLine=(ifs.readline()).rstrip()
                dummy = sLine.split(cDelimiter_in)
                ncolumn_out = len(dummy)
                #check ncolumn_in count
                if ncolumn_out < 1 :
                    print( 'The file has no data!' )  
                    pass

                aData_out = np.full( (nrow_out, ncolumn_out), '  ' , dtype=object )         
                aData_out[0] = dummy
                for iRow in range(1, nrow_out):
                    sLine=(ifs.readline()).rstrip()
                    if iFlag_iFlag_remove_quota ==1:
                        sLine.replace('"','')
                    else:
                        pass    

                    aDummy = sLine.split(cDelimiter_in)
                    iDummy = len(aDummy)
                    if iDummy == ncolumn_out:
                        aData_out[iRow] = aDummy
                    else:
                        #something is missing
                        aDummy2  =np.full(ncolumn_out, '  ' , dtype=object)
                        aDummy2[0: iDummy]= aDummy
                        aData_out[iRow]= aDummy2
                        pass
                
        else :
            if iFlag_column == 1:
                #check ncolumn_in count
                if ncolumn_out < 1 :                
                    return aData_out                      

                aData_out = np.full( (nrow_out, ncolumn_out), '  ', dtype=object )
                for iRow in range(nrow_out):
                    dummy1 = ifs.readline()
  
                    dummy2 = dummy1.rstrip()
                    if iFlag_iFlag_remove_quota ==1:
                        dummy2=dummy2.replace('"','')
                    else:
                        pass

                    dummy3 = dummy2.split()
        
                    aData_out[iRow] = dummy3
            else :
                sLine=(ifs.readline()).rstrip()
                dummy = sLine.split()
                ncolumn_out = len(dummy)
                if ncolumn_out < 1 :                
                    return aData_out          

                aData_out = np.full( (nrow_out, ncolumn_out), '  ', dtype=object )
                if iFlag_iFlag_remove_quota ==1:
                    sLine=sLine.replace('"','')
                else:
                    pass

                aData_out[0] =  dummy
                for iRow in range(1,nrow_out):
                    dummy1 = ifs.readline()
                    #print('New line is: ' + dummy1)
                    dummy2 = dummy1.rstrip()
                    if iFlag_iFlag_remove_quota ==1:
                        dummy2=dummy2.replace('"','')
                    else:
                        pass

                    dummy3 = dummy2.split()
               
                    aData_out[iRow] = dummy3


        ifs.close()

    else :
        print('file does not exist')

    return aData_out
