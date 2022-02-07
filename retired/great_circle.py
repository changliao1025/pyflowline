d1 = pVertex_start.calculate_distance(pVertex_in)
            d2 = pVertex_end.calculate_distance(pVertex_in)
            d3 = d1 +d2 -self.dLength

            angle3deg = calculate_angle_betwen_vertex(\
                 pVertex_start.dLongitude_degree, pVertex_start.dLatitude_degree,\
                 pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,\
                 pVertex_end.dLongitude_degree,pVertex_end.dLatitude_degree)



            if  angle3deg > 179 and d3 < 0.11: #care
                iFlag = 1
                dDistance = d1
            else:
                iFlag = 0