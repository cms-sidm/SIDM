from coffea.nanoevents import NanoAODSchema

class LLPNanoAODSchema(NanoAODSchema):
    """LLPNano schema builder
    LLPNano is an extended NanoAOD format that includes DSA Muons and improved displacement info
    """
    mixins = {
        **NanoAODSchema.mixins,
        "Muon": "LLPMuon", #Adds the matched_dsa_muon property on top of the normal NanoAOD Muon behavior
        "DSAMuon": "DSAMuon",
    }

    all_cross_references = {
        **NanoAODSchema.all_cross_references,
        "Muon_dsaMatch1idx": "DSAMuon",
        "Muon_dsaMatch2idx": "DSAMuon",
        "Muon_dsaMatch3idx": "DSAMuon",
        "Muon_dsaMatch4idx": "DSAMuon",
        "Muon_dsaMatch5idx": "DSAMuon",
        "DSAMuon_muonMatch1idx": "Muon",
        "DSAMuon_muonMatch2idx": "Muon",
        "DSAMuon_muonMatch3idx": "Muon",
        "DSAMuon_muonMatch4idx": "Muon",
        "DSAMuon_muonMatch5idx": "Muon",
    }
    nested_items = {
        **NanoAODSchema.nested_items,
        "Muon_dsaIdxG": [
            "Muon_dsaMatch1idxG",
            "Muon_dsaMatch2idxG",
            "Muon_dsaMatch3idxG",
            "Muon_dsaMatch4idxG",
            "Muon_dsaMatch5idxG",
        ],
        "DSAMuon_muonIdxG": [
            "DSAMuon_muonMatch1idxG",
            "DSAMuon_muonMatch2idxG",
            "DSAMuon_muonMatch3idxG",
            "DSAMuon_muonMatch4idxG",
            "DSAMuon_muonMatch5idxG",
        ],
    }

    @classmethod
    def behavior(cls):
        """Behaviors necessary to implement this schema (dict)"""

        import awkward
        import numpy
        from coffea.nanoevents.methods import nanoaod, base, candidate, vector
        from dask_awkward import dask_property


        # nanoaod.behavior.update(awkward._util.copy_behaviors("PtEtaPhiMCandidate", "DSAMuon", nanoaod.behavior))
        nanoaod.behavior.update(awkward._util.copy_behaviors("Muon", "DSAMuon", nanoaod.behavior))

        @awkward.mixin_class(nanoaod.behavior)
        class DSAMuon(candidate.PtEtaPhiMCandidate, base.NanoCollection, base.Systematic):
            """LLPNanoAOD DSA muon object"""
            @dask_property
            def matched_muons(self):

                """The matched PF muons (up to 5) as determined by the NanoAOD branch muonMatchNidx)"""
                muon_match_total = awkward.concatenate([
                self.muonMatch1[:, :, numpy.newaxis],
                self.muonMatch2[:, :, numpy.newaxis],
                self.muonMatch3[:, :, numpy.newaxis],
                self.muonMatch4[:, :, numpy.newaxis],
                self.muonMatch5[:, :, numpy.newaxis],
                ], axis=2) # Result: (events, dsa_muons, 5)

                pf_matches = self._events().Muon._apply_global_index(self.muonIdxG)

                concat = awkward.with_field(pf_matches, muon_match_total, where="numMatch")

                return concat


            @matched_muons.dask
            def matched_muons(self, dask_array):

                muon_match_total = awkward.concatenate([
                dask_array.muonMatch1[:, :, numpy.newaxis],
                dask_array.muonMatch2[:, :, numpy.newaxis],
                dask_array.muonMatch3[:, :, numpy.newaxis],
                dask_array.muonMatch4[:, :, numpy.newaxis],
                dask_array.muonMatch5[:, :, numpy.newaxis],
                ], axis=2) # Result: (events, dsa_muons, 5)

                pf_matches = dask_array._events().Muon._apply_global_index(dask_array.muonIdxG)

                concat = awkward.with_field(pf_matches, muon_match_total, where="numMatch")

                return concat

            # @property
            # def mass(self):
            #     return awkward.ones_like(self.pt)*0.106


        nanoaod._set_repr_name("DSAMuon")

        DSAMuonArray.ProjectionClass2D = vector.TwoVectorArray  # noqa: F821  # pylint: disable=undefined-variable
        DSAMuonArray.ProjectionClass3D = vector.ThreeVectorArray  # noqa: F821  # pylint: disable=undefined-variable
        DSAMuonArray.ProjectionClass4D = DSAMuonArray  # noqa: F821  # pylint: disable=undefined-variable
        DSAMuonArray.MomentumClass = vector.LorentzVectorArray  # noqa: F821  # pylint: disable=undefined-variable

        ##### Now update standard Muon class
        nanoaod.behavior.update(awkward._util.copy_behaviors("Muon", "LLPMuon", nanoaod.behavior))

        @awkward.mixin_class(nanoaod.behavior)
        class LLPMuon(nanoaod.Muon, base.NanoCollection, base.Systematic):
            """LLPNanoAOD Muon object"""
            @dask_property
            def matched_dsa_muons(self):

                """The matched PF muons (up to 5) as determined by the NanoAOD branch muonMatchNidx)"""
                muon_match_total = awkward.concatenate([
                self.dsaMatch1[:, :, numpy.newaxis],
                self.dsaMatch2[:, :, numpy.newaxis],
                self.dsaMatch3[:, :, numpy.newaxis],
                self.dsaMatch4[:, :, numpy.newaxis],
                self.dsaMatch5[:, :, numpy.newaxis],
                ], axis=2) # Result: (events, dsa_muons, 5)

                dsa_matches = self._events().DSAMuon._apply_global_index(self.dsaIdxG)

                concat = awkward.with_field(dsa_matches, muon_match_total, where="numMatch")

                return concat


            @matched_dsa_muons.dask
            def matched_dsa_muons(self, dask_array):

                muon_match_total = awkward.concatenate([
                dask_array.dsaMatch1[:, :, numpy.newaxis],
                dask_array.dsaMatch2[:, :, numpy.newaxis],
                dask_array.dsaMatch3[:, :, numpy.newaxis],
                dask_array.dsaMatch4[:, :, numpy.newaxis],
                dask_array.dsaMatch5[:, :, numpy.newaxis],
                ], axis=2) # Result: (events, dsa_muons, 5)

                dsa_matches = dask_array._events().DSAMuon._apply_global_index(dask_array.dsaIdxG)

                concat = awkward.with_field(dsa_matches, muon_match_total, where="numMatch")

                return concat

        nanoaod._set_repr_name("LLPMuon")

        LLPMuonArray.ProjectionClass2D = vector.TwoVectorArray  # noqa: F821  # pylint: disable=undefined-variable
        LLPMuonArray.ProjectionClass3D = vector.ThreeVectorArray  # noqa: F821  # pylint: disable=undefined-variable
        LLPMuonArray.ProjectionClass4D = LLPMuonArray  # noqa: F821  # pylint: disable=undefined-variable
        LLPMuonArray.MomentumClass = vector.LorentzVectorArray  # noqa: F821  # pylint: disable=undefined-variable



        return nanoaod.behavior
