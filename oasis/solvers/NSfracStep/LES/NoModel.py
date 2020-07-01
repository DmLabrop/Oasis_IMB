from dolfin import Constant

__all__ = ['les_setup', 'les_update']

def les_setup(**NS_namespace):
    return dict(nut_=Constant(0))


def les_update(**NS_namespace):
    pass
