def calculate_potentiometric(dElevation_bedrock_in, dThickness_ice_in, dDensity_ice_in=None):
    #calculate potentiometric surface

    #if both are missing value, keep it as missing
    if dElevation_bedrock_in == -9999 or dThickness_ice_in == -9999:
        return -9999

    if dDensity_ice_in is None:
        dDensity_ice = 917.0 #unit
    else:
        dDensity_ice = dDensity_ice_in

    dDensity_water = 1000.0 #unit

    iOption = 1
    dGravity = 9.81
    if iOption == 1:
        #reference
        #https://www.tandfonline.com/doi/full/10.1080/15481603.2016.1230084
        dRatio_k = 1.0
        dDummy0 = dDensity_ice / dDensity_water
        dDummy1 = dRatio_k * dDummy0 * dThickness_ice_in
        dPotentiometric = dDummy1 + dElevation_bedrock_in

    else:
        #unit and quality check?
        #reference
        #use equation 1 from https://onlinelibrary.wiley.com/doi/epdf/10.1002/hyp.7343
        dFactor = 0.1
        dElevation_ice =  dElevation_bedrock_in + dThickness_ice_in
        dDummy=  dElevation_ice + dElevation_bedrock_in * dFactor
        dPotentiometric = dDensity_ice * dGravity *  dDummy

    return dPotentiometric


if __name__ == "__main__":

    dElevation_bedrock = 100.0
    dThickness_ice = 1000.0
    dDensity_ice = 917.0

    dPotentiometric = calculate_potentiometric(dElevation_bedrock, dThickness_ice, dDensity_ice)

    print(dPotentiometric)