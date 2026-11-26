"""Module to define the Histogram and Axis classes"""

import traceback
import hist
import awkward as ak


class Histogram:
    """Class to represent histograms

    Histogram mostly exists so that hist.Hists and the appropriate filling arguments can be
    defined in one place. In addition to the filling function associated with each Axis, the user
    can optionally provide an event mask that is applied to all object collections used to fill
    the histogram, e.g. to ensure only events with >=2 muons are used to fill dR(mu, mu) hists.
    """

    def __init__(self, axes, storage="weight", evt_mask=None, entry_factors=None):
        """entry_factors, if given, is (categories, function): function(objs, mask) returns
        {category: factor per entry}, jagged like the values of the last axis. The hist then has
        an "sf" axis with these categories, and each entry is filled once per category with its
        event weight times its factor. This is for weights that belong to the objects a hist
        counts rather than to the event, such as the scale factors of the lepton jet that a probe
        of a control sample is matched to."""
        self.axes = axes
        self.storage = storage
        # Allow all events to pass if no mask is explicitly provided
        self.evt_mask = (lambda objs: slice(None)) if evt_mask is None else evt_mask
        self.entry_factors = entry_factors
        self.hist = None
        self.systematics = None

    @classmethod
    def simple_hist(cls, obj, attr, absval, nbins, xmin, xmax, label):
        """Method to simplify creation of basic obj.attribute hists"""
        # define fill function
        def f(objs, mask):
            if attr == "n":
                return ak.num(objs[obj])
            elif absval:
                return abs(objs[obj].__getattr__(attr))
            else:
                return objs[obj].__getattr__(attr)
        return cls([
            Axis(hist.axis.Regular(nbins, xmin, xmax, name=f"{obj}_{attr}", label=label), f)
            ])

    def make_hist(self, name, channels=None, lj_reco_choices=None, systematics=None):
        """Build associated hist.Hist

        Perform outside __init__ because channels aren't known until runtime.
        systematics optionally adds a "syst" axis with one category per weight
        variation, "nominal" included.
        """
        self.name = name

        # optionally add channels axis to hist
        if channels is not None:
            channel_axis = hist.axis.StrCategory(channels, name="channel")
            self.axes = [Axis(channel_axis, lambda objs, mask: objs["ch"])] + self.axes

        # optionally add lj_reco axis to hist
        if lj_reco_choices is not None:
            lj_reco_axis = hist.axis.StrCategory(lj_reco_choices, name="lj_reco")
            self.axes = [Axis(lj_reco_axis, lambda objs, mask: objs["lj_reco"])] + self.axes

        axes = [a.axis for a in self.axes]

        # optionally add a syst axis after the category axes. It has no fill
        # function: fill() sets it once per variation.
        n_category_axes = (channels is not None) + (lj_reco_choices is not None)
        if systematics is not None:
            self.systematics = list(systematics)
            axes.insert(n_category_axes, hist.axis.StrCategory(self.systematics, name="syst"))
            n_category_axes += 1

        # the sf axis of entry_factors, which fill() also sets once per category
        if self.entry_factors is not None:
            axes.insert(n_category_axes, hist.axis.StrCategory(list(self.entry_factors[0]), name="sf"))

        self.hist = hist.Hist(*axes, storage=self.storage)

    def fill(self, objs, evt_weights, verbose=False):
        """Fill associated hist.Hist

        evt_weights is one weight per event or {syst category: weights}. A hist
        with a syst axis is filled once per category given; categories left out
        stay empty, which is how data fills only "nominal". A hist without one
        uses the "nominal" weights.
        """
        # Create fill args, warning user and skipping hists that cannot be filled
        try:
            mask = self.evt_mask(objs)
            fill_args = {a.name: a.fill_func(objs, mask) for a in self.axes}
            factors = ({None: None} if self.entry_factors is None
                       else self.entry_factors[1](objs, mask))
        except (AttributeError, KeyError, ValueError) as e:
            print(f"Warning: the fill function of a histogram with the name {self.name} could not be evaluated and will be skipped ({e!r})")
            if verbose:
                print(traceback.format_exc())
            return

        if not isinstance(evt_weights, dict):
            evt_weights = {"nominal": evt_weights}
        if self.systematics is None:
            weight_sets = {None: evt_weights["nominal"]}
        else:
            weight_sets = {s: evt_weights[s] for s in self.systematics if s in evt_weights}

        # Use last axis to define weight structure to avoid channels axis
        shape = fill_args[self.axes[-1].name]
        for name in fill_args.keys():
            if name not in ("channel", "lj_reco"):
                fill_args[name] = ak.flatten(fill_args[name], axis=None)

        for syst, weights in weight_sets.items():
            if not self._fill_categories(fill_args, weights[mask]*ak.ones_like(shape), factors, syst):
                return

    def _fill_categories(self, fill_args, weights, factors, syst):
        """Fill the hist once per category of factors (once, without an sf axis), each entry
        weighted by its weight times its factor; False if it cannot be filled"""
        for category, factor in factors.items():
            args = dict(fill_args)
            args["weight"] = ak.flatten(weights if factor is None else weights*factor, axis=None)
            if syst is not None:
                args["syst"] = syst
            if category is not None:
                args["sf"] = category
            # Fill hist, warning user and skipping hists that cannot be filled
            try:
                self.hist.fill(**args)
            except ValueError:
                print(f"Warning: a histogram with the name {self.name} could not be filled and will be skipped")
                return False
        return True

class Axis:
    """Class to represent histogram axes

    Axis just bundles together hist.axis objects and functions to fill them.
    """

    def __init__(self, axis, fill_func):
        self.axis = axis
        self.name = self.axis.name
        self.fill_func = fill_func
