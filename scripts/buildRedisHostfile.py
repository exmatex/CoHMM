import os;

def genNodeList():
    #Get environmental variable
    #nodeString = os.environ['SLURM_JOB_NODELIST']
    #For debug/devel
    nodeString = "cn[3,46,52,61,111-113,116-118]"
    #Remove leading 'cn'
    nodeString = nodeString.lstrip('cn')
    nodeList = []
    if (nodeString.startswith('[') == False):
        #Just the one
        nodeList.append(int(nodeString))
    else:
        #Strip the brackets, we no longer need them
        nodeString = nodeString.strip('[]')
        #Tokenize the string
        tokenList = nodeString.rsplit(',')
        #Iterate over tokens
        for token in tokenList:
            #Two cases: Single or range
            if '-' in token:
                inclusiveRange = token.rsplit('-')
                start = int(inclusiveRange[0])
                end = int(inclusiveRange[1]) + 1
                for node in xrange(start, end):
                    nodeList.append(node)
            else:
                nodeList.append(int(token))
    return nodeList;

def assignServers(nodeList, numServers):
    serverList = []
    #See if this is one of the easy case
    if (numServers == 1):
        #Just one server, so return the first
        serverList.append(nodeList[0])
    elif (numServers >= len(nodeList)):
        #Everyone is a server
        serverList = nodeList
    else:
        #We need to assign servers
        #TODO
        serverList = []
    return serverList


def main():
    nodeList = genNodeList();
    serverList = assignServers(nodeList, 4)
    print nodeList

if __name__ == "__main__":
    main()
