                if iRow >= 1:#0
                    lCellID0 = lCellID_center - 1
                    aNeighbor.append(lCellID0)

                if iColumn < ncolumn :#1
                    if iColumn %2 ==0:
                        lCellID1 = nrow * iColumn + iRow 
                        aNeighbor.append(lCellID1)
                    else:
                        if iRow != 1:
                            lCellID1 = nrow * iColumn + iRow 
                            aNeighbor.append(lCellID1)

                if iColumn < ncolumn :#2
                    if iColumn %2 ==1:
                        lCellID2 = nrow * iColumn + iRow
                        aNeighbor.append(lCellID2)
                    else:
                        if iRow!=nrow:
                            lCellID2 = nrow * iColumn + iRow
                            aNeighbor.append(lCellID2)
                        
                if iRow < nrow:#3
                    lCellID3 = lCellID_center + 1
                    aNeighbor.append(lCellID3)

                if iColumn > 1 : #4
                    if iColumn %2 ==1:
                        lCellID4 = nrow * iColumn + iRow 
                        aNeighbor.append(lCellID4)
                    else:
                        if iRow!=nrow:
                            lCellID4 = nrow * iColumn + iRow 
                            aNeighbor.append(lCellID4)


                if iColumn > 1 :#5
                    if iColumn %2 ==0:
                        lCellID5 = nrow * iColumn + iRow 
                        aNeighbor.append(lCellID5)
                    else:
                        if iRow !=1:
                            lCellID5 = nrow * iColumn + iRow 
                            aNeighbor.append(lCellID5)