#!/usr/bin/env python3
"""
Wrapper to the PHYEX sei
"""

# Copyright (c) Météo France (2025-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

# pylint: disable=invalid-name

# On importe les librairies nécessaires
import os
import tempfile
import numpy
import f90nml
import numpy as np
from tqdm import tqdm

# On importe ici les classes extèrieures
from condi_init import InitialCond


from pyphyex import PYICE4_SEDIMENTATION, PYLIMA_SEDIMENTATION, PYINI_PHYEX

class Phyex():
    """
    Generic class to call the different sedimentation schemes of PHYEX
    """
    def __init__(self, method, number_stitches, number_bin,
                 delta_t, speed_max, esp, CFL,type_init):
        """
        Here we initialise the non-spatial fixed parameters and allow important variables 
        to travel between functions. We also call the initialisation.
        """
        self.method = method
        self.number_stitches = number_stitches
        self.length_sim = 5000  # length of simulation in seconds
        self.delta_t = delta_t # length of time step in seconds
        self.nb_step = self.length_sim // self.delta_t  # number of time step
        #self.nb_diam = number_bin # number of type of diameter
        self.esp = esp

        # Geometry and initial content
        condi_init = InitialCond(self.number_stitches, esp, "bulk", nb_classes=1,mode=type_init)
        self.N_profile = np.array(condi_init.data["concentration"])
        self.rho_r_profile = np.array(condi_init.data["rho_r"])

        self.vertical_boundaries = condi_init.levels_boundaries
        self.levels = condi_init.grid

        # Initialisation of variables
        self.water_on_floor = 0.  # Mass of water wich has touched the ground
        self.wat_flo_on_time=[0.]  # List of this mass in time, here, time = 0 soit no precipitation

        # PHYEX init
        if method == 'EULE':
            self.ccloud = 'ICE3'
            nml = {'NAM_PARAM_ICEn': {'LSEDIC':True,
                                      'CSEDIM':'SPLI',
                                      'XSPLIT_MAXCFL':.8}}
        elif method == 'EULE2':
            self.ccloud = 'LIMA'
            nml = {'NAM_PARAM_LIMA': {'LSEDC':True,
                                      'LSEDI':True,
                                      'NMOM_C':2,
                                      'NMOM_R':2,
                                      'NMOM_I':2,
                                      'NMOM_S':2,
                                      'NMOM_G':2,
                                      'NMOM_H':2},
                   'NAM_NEBn': {'LCONDBORN': True}}
        elif method == 'STAT':
            self.ccloud = 'ICE3'
            nml = {'NAM_PARAM_ICEn': {'LSEDIC':True,
                                      'CSEDIM':'STAT'}}
        with tempfile.NamedTemporaryFile() as namel:
            self.nml = f90nml.read(namel.name)
            for k, v in nml.items():
                self.nml[k] = v
            self.nml.write(namel.name, force=True)
            PYINI_PHYEX('AROME', 33, namel.name, False, 20, 0, 1, self.delta_t, 20.,
                        self.ccloud, 'EDKF', 'TKEL',
                        LDCHANGEMODEL=True, LDDEFAULTVAL=True, LDREADNAM=True, LDCHECK=True,
                        KPRINT=0, LDINIT=True)
            os.remove('fort.20')

    def run(self):
        """
        On fait tourner le modèle pour chaque pas de temps.
        """
        NIJT = 1
        NKT = len(self.levels)
        PRHODREF = numpy.ones((NKT, NIJT))
        P = numpy.ones((NKT, NIJT)) * 101315.
        Theta = numpy.ones((NKT, NIJT)) * 300.
        T = numpy.ones((NKT, NIJT)) * 300.
        boundaries = numpy.array([0] + list(self.vertical_boundaries))
        dzz = (boundaries[1:] - boundaries[:-1]).reshape((NKT, 1))
        PLVFACT = PLSFACT = numpy.ones((NKT, NIJT))
        rhodj = dzz * PRHODREF
        PTHS = Theta / self.delta_t
        KRR = 6
        PRT = numpy.ndarray((KRR, NKT, NIJT))
        PCT = numpy.ndarray((KRR, NKT, NIJT))
        for k in range(KRR):
            PRT[k] = self.rho_r_profile.reshape((NKT, 1))
            PCT[k] = self.N_profile.reshape(NKT, 1)
        PRS = PRT / self.delta_t
        PCS = PCT / self.delta_t
        PQCT = PQRT = PQIT = PQST = PQGT = numpy.zeros((NKT, NIJT))
        PQCS = PQRS = PQIS = PQSS = PQGS = PEFIELDW = numpy.zeros((NKT, NIJT))
        PCONC3D = numpy.ones((NKT, NIJT)) * 3.E8

        rho_r_profile = []
        ct_profile = []

        print (" ")
        print ("---------------------------------------")

        format_affichage_time = "{l_bar}|{bar}| {n_fmt}/{total_fmt} pas de temps | temps écoulé : {elapsed} < temps restant {remaining} |" # definition of the format of the tqdm bar
        format_affichage_diam = "-> {desc} |{bar}| {n_fmt}/{total_fmt} bins  | {elapsed}<{remaining} |                                                                                           "

        cumul = [0. for spec in ['c', 'r', 'i', 's', 'g']]

        #print(f'content before (max r is {self.rho_r_profile.max()})')
        #for i, spec in enumerate(['v', 'c', 'r', 'i', 's', 'g']):
        #    print(spec, (PRT[i]*dzz*PRHODREF).sum())
        #print()
        for _ in tqdm(range(self.nb_step), bar_format = format_affichage_time, desc = f"Avancement total {self.ccloud} ", position = 0, colour = 'blue'):
            inst = numpy.zeros((KRR, NIJT))
            if self.ccloud == 'LIMA':
                for i, spec in enumerate(['c', 'r', 'i', 's', 'g']):
                    # PCT/PRT/PCS/PRS are declared with KRR (first index is for vapor)
                    # cloud (i=0) is the second array, indexed 1 in python and 2 in FORTRAN
                    KISHAPE = 1
                    result = PYLIMA_SEDIMENTATION(NIJT, NKT, 1, 0,
                                                  'L' if spec in ('c', 'r') else 'I',
                                                  self.nml['NAM_PARAM_LIMA'][f'NMOM_{spec.upper()}'],
                                                  i + 2, KISHAPE, 0, self.delta_t, False, dzz,
                                                  PRHODREF, 0., P, T, numpy.zeros((NKT, NIJT)),
                                                  numpy.zeros((NKT, NIJT)), PRS[i + 1], PCS[i+ 1],
                                                  missingOUT=['PQS'])
                    (_, _, _, PRS[i + 1], PCS[i + 1], inst[i + 1], _) = result
            else:
                result = PYICE4_SEDIMENTATION(NIJT, NKT, 1, 0, False, False, self.delta_t, KRR,
                                              dzz, 0., PLVFACT, PLSFACT,
                                              PRHODREF, P, Theta,
                                              T, rhodj, PTHS, PRT, PRS, PCONC3D,
                                              PQCT, PQRT, PQIT, PQST, PQGT,
                                              PQCS, PQRS, PQIS, PQSS, PQGS,
                                              PEFIELDW, 0, PSEA=numpy.zeros((NIJT,)),
                                              PTOWN=numpy.zeros((NIJT,)),
                                              missingOUT=['PFPR', 'PINPRH'])
                (PTHS, PRS, inst[1], inst[2], inst[3],
                 inst[4], _, _, _, _, _) = result
            PCT = PCS * self.delta_t
            PRT = PRS * self.delta_t
            for i, spec in enumerate(['c', 'r', 'i', 's', 'g']):
                cumul[i] += inst[i+1, 0] * self.delta_t * 1000.
                if self.esp == spec:
                    rho_r_profile.append(PRT[i+1, :, 0].copy())
                    ct_profile.append(PCT[i+1, :, 0].copy())
                    self.wat_flo_on_time.append(inst[i+1, 0] * self.delta_t * 1000.)  # m/s -> kg/m2
        #print('Content after (in atmosphere, at surface, total):')
        #for i, spec in enumerate(['v', 'c', 'r', 'i', 's', 'g']):
        #    print(spec, (PRT[i]*dzz*PRHODREF).sum(), cumul[i-1] if i !=0 else 0, (PRT[i]*dzz*PRHODREF).sum() + (cumul[i-1] if i !=0 else 0))

        self.rho_r_profile = rho_r_profile
        self.ct_profile = ct_profile

        print ("---------------------------------------")

        return self.wat_flo_on_time, self.rho_r_profile, self.ct_profile


class Eule(Phyex):
    """
    Eulerian 1-moment scheme
    """
    def __init__(self, *args, **kwargs):
        super().__init__('EULE', *args, **kwargs)

class Eule2(Phyex):
    """
    Eulerian 2-moment scheme
    """
    def __init__(self, *args, **kwargs):
        super().__init__('EULE2', *args, **kwargs)

class Stat(Phyex):
    """
    Statistic scheme
    """
    def __init__(self, *args, **kwargs):
        super().__init__('STAT', *args, **kwargs)
