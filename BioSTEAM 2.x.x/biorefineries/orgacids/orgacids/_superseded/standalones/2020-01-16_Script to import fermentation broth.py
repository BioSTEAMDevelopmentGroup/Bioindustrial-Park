# The orgacids module is more or less completed aside from separation processes
# You can check flow info, system setup, export diagrams/reports, etc.

# To import fermentation broth from orgacids module
import biosteam as bst
from orgacids import system

broth_for_separation = bst.find.stream.broth_for_separation