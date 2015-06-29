import os;

def genRedisServers(numServers):
    numNodes = os.environ('FOO')
    numNodes = int(numNodes)
    nodeString = os.environ('MAN')
    #cn[3,46,52,61,111-113,116-118] is an example
    #Build nodeList from nodeString
