from .data_structures import Node
from .data_structures import TreeWalker
from slugify import slugify


class GPRObject(Node):

    def __init__(self, name, *args):
        Node.__init__(self, name)
        self.args = args

    def to_command(self):
        s = self.fs.format(*self.args)
        return s

    def __str__(self):
        return self.to_command()


class GPRObjectCreator:

    def __init__(self):
        self.types = {
            'discretisation': '#dx_dy_dz: {} {} {}',
            'time_window': '#time_window: {}',
            'title': '#title: {}',
            'edge': '#edge: {} {} {} {} {} {} {}',
            'box': '#box: {} {} {} {} {} {} {}',
            'domain': '#domain: {} {} {}',
            'waveform': '#waveform: {} {} {} {}',
            'voltage_source': '#voltage_source: {} {} {} {} {} {}',
            'cylinder': '#cylinder: {} {} {} {} {} {} {} {}',
            'hertzian_dipole': '#hertzian_dipole: {} {} {} {} {} {} {}',
            'snapshot': '#snapshot: {} {} {} {} {} {} {} {} {} {} snapshot{}',
            'rx': '#rx: {} {} {}',
            'geometry_view': '#geometry_view: {} {} {} {} {} {} {} {} {} {} {}',
            'dx_dy_dz': '#dx_dy_dz: {} {} {}',
            'material': '#material: {} {} {} {} {}',
            'time_window': '#time_window: {}',
            'wrapper': 'wrapper',
            'waveform': '#waveform: {} {} {} {}',
            'transmission_line': '#transmission_line: {} {} {} {} {} {}'
        }

    def create(self, name, *args):
        fs = self.types.get(name, None)
        if fs is None:
            raise Exception('Unknown GPRObject Type: ', name)
        if fs.count('{}') != len(args):
            raise Exception('Incorrect number of arguments to create: ', name)
        e = GPRObject(name, *args)
        e.fs = fs
        return e


class Scene(Node):

    def __init__(self):
        Node.__init__(self, 'scene')

    def to_commands(self):
        s = ''
        tw = TreeWalker()
        nodes = tw.getBreadthFirstNodes(self)
        for node in nodes:
            if node.name != 'wrapper':
                s += node.to_command() + '\n'
        return s

    def getTitle(self):
        tw = TreeWalker()
        t = tw.getNode(self, 'title').args[0]
        return t

    def write(self):
        filepath = self.getTitle() + '.in'
        f = open(filepath, 'w')
        f.write(self.to_commands())
        f.close()
        return filepath


def write_scene(scene):
    commands = scene.to_commands()
    t = scene.getTitle()
    fp = slugify(t) + '.in'
    with open(fp, 'w') as f:
        f.write(commands)
    return fp
