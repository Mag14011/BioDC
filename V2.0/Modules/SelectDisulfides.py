################################################################################################################################################

def SelectDisulfides(InputDict, LaunchDir):
    while True:
        try:
            if ("NumDisulfide" in InputDict):
                NumDisulfide = int(InputDict["NumDisulfide"])
                break
            else:
                NumDisulfide = int(input(" How many disulfide linkages are present? "))
                break
        except ValueError:
            print(" Your entry must be an integer.")

        print(f"NumDisulfide = {NumDisulfide}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

    if (NumDisulfide != 0):
        print(" For each disulfide bond, please enter the pair of Cys residue IDs.")
        if ("SelPairIDs" in InputDict):
            DisulfResList = InputDict["DisulfResList"]
            DisulfPairID = InputDict["DisulfPairID"]
        else:
            idx = 0
            DisulfResList = " "
            DisulfPairID = [0]*NumDisulfide
            for n in range(NumDisulfide):
                SelPairIDs = input(f"  Disulfide-linked Cys pair {idx+1}: ")

                print(SelPairIDs, file=open('DisulfideDefinitions.txt', 'w'))
                DisulfResList += SelPairIDs+" "
                DisulfPairID[idx] = list(map(int,SelPairIDs.split()))

                idx+=1
        print(f"SelPairIDs = {DisulfResList}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

    return DisulfResList, DisulfPairID

################################################################################################################################################
