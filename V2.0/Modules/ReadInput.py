################################################################################################################################################
#import iteritems

################################################################################################################################################
def ReadInput():

    InputDict = {}
    with open("input.txt") as inp:
        Lines_inp = inp.readlines()
        for line in Lines_inp:
            EntryLength = len(line.strip().split(" "))

            #For keywords that have one and only one value
            if (EntryLength <= 3):
                keyword = line.strip().split(" ")[0]

                try:
                    InputDict[keyword] = line.strip().split(" ")[2]
                except IndexError:
                    InputDict[keyword] = ''

                #If keyword = MutResArray for one mutaiton
                if (keyword == "MutResArray"):
                    MutID = 0
                    MutResArray = [0]*(EntryLength-2)

                    for i in range(2, EntryLength):
                        MutResArray[MutID] = line.strip().split(" ")[i].split("-")
                        MutID+=1
                    InputDict["MutResArray"] = MutResArray

            #For keywords with more than one value
            if (EntryLength > 3):
                keyword = line.strip().split(" ")[0]

                value=" "
                for i in range(2, (EntryLength)):
                    value = f"{value} "+line.strip().split(" ")[i]
                InputDict[keyword] = value

                #If multiple mutations are selected
                if (keyword == "MutResArray"):
                    MutID = 0
                    MutResArray = [0]*(EntryLength-2)

                    for i in range(2, EntryLength):
                        MutResArray[MutID] = line.strip().split(" ")[i].split("-")
                        MutID+=1
                    InputDict["MutResArray"] = MutResArray

                #If disulfide linkages are indicated
                if (keyword == "SelPairIDs"):
                    try:
                        PairNum = int((EntryLength-2)/2)
                        DisulfPairID = [0]*(PairNum)
                    except TypeError:
                        print("PairNum is not an integer.")

                    PairNumID = 0
                    DisulfResList = ''
                    for i in range(2, (EntryLength-1)):
                        if (i % 2 == 0):
                            SelPairIDs = f"{line.strip().split(' ')[i]} {line.strip().split(' ')[i+1]}"
                            DisulfResList += SelPairIDs+" "
                            DisulfPairID[PairNumID] = list(map(str,SelPairIDs.split()))
                            PairNumID+=1
                
                    InputDict["DisulfResList"] = DisulfResList
                    InputDict["DisulfPairID"] = DisulfPairID

                #If multiple redox states are specified because there are multiple hemes
                if (keyword == "RedoxState"):
                    idx = 0
                    RedoxState = [0]*(EntryLength-2)

                    for i in range(2, EntryLength):
                        RedoxState[idx] = line.strip().split(" ")[i]
                        idx+=1
                    InputDict["RedoxState"] = RedoxState

    return InputDict
################################################################################################################################################

#InputDict = ReadInput()
#for x in InputDict:
#    print(repr(x),":",InputDict[x])

