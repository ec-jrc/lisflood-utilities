import os


def pcraster_command(cmd, files=None):

    if files is not None:
        # replace placeholders in command with actual filenames
        for alias, realname in files.items():
            cmd = cmd.replace(alias, '"{}"'.format(realname))

    os.system(cmd)
    return cmd
