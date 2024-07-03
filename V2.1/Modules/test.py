
def ReadInput():

    InputDict = {}
    with open("input.txt") as inp:
        Lines_inp = inp.readlines()
        for line in Lines_inp:
            EntryLength = len(line.strip().split(" "))

            if (EntryLength <= 3):
                try:
                    InputDict[line.strip().split(" ")[0]] = line.strip().split(" ")[2]
                except IndexError:
                    InputDict[line.strip().split(" ")[0]] = ''

            if (EntryLength > 3):
                keyword = line.strip().split(" ")[0]

                value=" "
                for i in range(2, (EntryLength)):
                    value = f"{value} "+line.strip().split(" ")[i]
                InputDict[keyword] = value

    return InputDict




