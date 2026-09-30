"""Module to define the Selection and JaggedSelection classes"""

# columnar analysis
import awkward as ak
# local
from sidm.definitions.cuts import evt_cut_defs, obj_cut_defs
from sidm.tools.cutflow import SimpleCutflow
import traceback

class Selection:
    """Class to represent the collection of cuts that define a Selection

    A selection consists of event-level cuts which reject whole events.

    All available cuts are defined in sidm.definitions.cuts. The specific cuts that define each
    selection are accepted by Selection() as lists of strings.
    """

    def __init__(self, cuts, verbose=False):
        self.evt_cuts = cuts # list of names of cuts to be applied
        self.cutflow = SimpleCutflow() # will be filled later when cuts are evaluated
        self.verbose = verbose

    def apply_evt_cuts(self, objs):
        """Evaluate all event cuts and apply results to object collections"""
        self.cutflow.add_row("None", ak.count(objs["evt_weights"]), ak.sum(objs["evt_weights"]))
        sel_objs = objs.copy()
        for cut in self.evt_cuts:
            # evaluate cuts
            if self.verbose:
                print("Applying cut:", cut)
            try:
                mask = evt_cut_defs[cut](sel_objs)
            except Exception as e:
                # An I/O failure means an unreadable chunk, not an unevaluable cut.
                if isinstance(e, OSError) or "I/O operation" in str(e):
                    raise
                # A collection absent from this sample (a gen-level cut on data)
                # means the cut does not apply here, not that it is broken;
                # report it the way JaggedSelection.apply_obj_cuts already does
                # for object cuts on a missing collection.
                if (isinstance(e, KeyError) and e.args and isinstance(e.args[0], str)
                        and e.args[0] != cut and e.args[0] not in sel_objs):
                    print(f"Warning: {e.args[0]} not found in sample. "
                          f"The following cuts will not be applied: ['{cut}']")
                    continue
                print(f"Warning: Unable to evaluate {cut} Skipping.",e)
                # Skip this cut entirely. Event cuts were once collected and
                # applied together as a single AND after the loop, where falling
                # through here was harmless. They are now evaluated and applied
                # one at a time against the progressively slimmed collections,
                # so without this the loop below would re-apply the PREVIOUS
                # cut's mask and log a cutflow row under this cut's name.
                continue

            # apply to all object collections
            for name, obj in sel_objs.items():
                try:
                    sel_objs[name] = obj[mask]
                    if name == "evt_weights":
                        self.cutflow.add_row(cut, ak.count(sel_objs[name]), ak.sum(sel_objs[name]))
                except Exception as e:
                    if isinstance(e, OSError) or "I/O operation" in str(e):
                        raise
                    print(f"Warning: Unable to apply event cuts to {name}. Skipping. "
                          f"({type(e).__name__}: {e})")

        return sel_objs


class JaggedSelection:
    """Class to represent the collection of cuts that define a JaggedSelection

    A JaggedSelection consists of object-level cuts (for example, electron or lepton-jet-level cuts).
    Object-level cuts slim object collections and are stored as a dictionary of masks.

    All available cuts are defined in sidm.definitions.cuts. The specific cuts that define each
    selection are accepted by JaggedSelection() as lists of strings.
    """

    def __init__(self, cuts, verbose=False):
        self.obj_cuts = cuts # dict of cuts to be applied
        self.verbose = verbose

    def apply_obj_cuts(self, objs):
        """Apply object cuts sequentially"""
        sel_objs = objs.copy()
        for obj, cuts in self.obj_cuts.items():
            if obj not in objs:
                print(f"Warning: {obj} not found in sample. "
                      f"The following cuts will not be applied: {cuts}")
                continue

            for cut in cuts:
                if self.verbose:
                    print(f"Applying {obj} {cut}")
                try:
                    if obj == "muons" or obj == "dsaMuons":
                        # Notice the change in syntax here, driven by the need to apply the nested selection in the same way
                        # Eventually should change every cut to this form
                        sel_objs[obj] = sel_objs[obj][obj_cut_defs[obj][cut](sel_objs,sel_objs[obj])]
                    else:
                        sel_objs[obj] = sel_objs[obj][obj_cut_defs[obj][cut](sel_objs)]
                except Exception as e:
                    if isinstance(e, OSError) or "I/O operation" in str(e):
                        raise
                    print(f"Warning: Unable to apply {cut} for {obj}. Skipping. Error message: {e}")
                    traceback.print_exc()
        return sel_objs


class NestedSelection:
    """Class to represent the collection of cuts that define a NestedSelection

    A NestedSelection consists of object-level cuts to the matched objects
    (for example, the pf muons matched to dsa muons).

    All available cuts are defined in sidm.definitions.cuts. The specific cuts that define each
    selection are accepted by NestedSelection() as lists of strings.
    """

    def __init__(self, cuts, verbose=False):
        self.obj_cuts = cuts # dict of cuts to be applied
        self.verbose = verbose

    def apply_obj_cuts(self, all_objs, nested_objs, obj_name):
        """Apply object cuts sequentially"""
        # For example, all_objs would be the full sel_objs in the processor
        # and nested_objs could be sel_objs["dsaMuons"].matched_muons in the processor
        # Then obj_name would be "muons", because we want the *muon* cuts applied to these objects

        sel_nested_objs = nested_objs

        for obj, cuts in self.obj_cuts.items():
            if obj == obj_name:

                if self.verbose: print(f"Applying cuts to the nested {obj} collection")

                for cut in cuts:
                    try:
                        sel_nested_objs = sel_nested_objs[obj_cut_defs[obj][cut](all_objs, sel_nested_objs)]
                    except Exception as e:
                        print(f"Warning: Unable to apply {cut} for nested {obj_name} collection. Skipping.... {e}")
                return sel_nested_objs



        if self.verbose:
            print(f"Warning: Did not find cuts for {obj_name} in the config. "
                      f"No cuts will be applied to the nested objects")
        return sel_nested_objs

