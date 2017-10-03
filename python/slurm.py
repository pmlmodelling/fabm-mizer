import re

def getNodes(servers):
    if isinstance(servers, basestring):
        match = re.match(r'(.*)\[(.*)\](.*)', servers)
        if match is not None:
            # Hostnames in PBS/SLURM notation, e.g., node[01-06]
            servers = []
            left, middle, right = match.groups()
            for item in middle.split(','):
                if '-' in item:
                    start, stop = item.split('-')
                    for i in range(int(start), int(stop)+1):
                        servers.append('%s%s%s' % (left, str(i).zfill(len(start)), right))
                else:
                    servers.append('%s%s%s' % (left, item, right))
        else:
            # Comma-separated hostnames
            servers = servers.split(',')
        return tuple(servers)
    else:
        assert servers is None
        return ()
