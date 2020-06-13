def getBackwardName(forwardName):
    if (forwardName[-6:] != "fw.npy"):
        print("invalid input name")
        exit(0)
    backwardName = forwardName[:-6] + "bw.npy"
    return backwardName