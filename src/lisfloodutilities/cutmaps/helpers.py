import os


def pcraster_command(cmd='', files=None):

    if files is not None:
        # replace placeholders in command with actual filenames
        aliases = []
        realnames = []
        for alias, realname in files.items():
            aliases.append(alias)
            realnames.append(realname)

        for i, alias in enumerate(aliases):
            realname = realnames[i]
            cmd = cmd.replace(alias, '"{}"'.format(realname))

    os.system(cmd)
    return cmd
