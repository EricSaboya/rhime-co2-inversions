# *****************************************************************************
# nuclear14c.py
# Created: 30 Jan. 2025
# Author: Eric Saboya, School of Geographical Sciences, University of Bristol
# *****************************************************************************
# Radiocarbon nuclear power plant functions 
# *****************************************************************************

import os
import openghg
import numpy as np

def get_npp_eu_radd(start_date: str | None = None,
                    end_date: str | None = None, 
                   ):
    """
    Retrieve European Commission RAdioactive Discharges Database (RADD) 
    14CO2 reported annual values. 

    Does not include UK nuclear power plants.
    Args:
        start_date (str):
            Date from which to gather data
        end_date (str):
            Date until which to gather data

    Returns
    NPP 14C flux field
    """
    npp_14c = openghg.retrieve.get_flux(store="radiocarbon_store",
                                        species="dco2c14",
                                        source="radd-eu",
                                        domain="EUROPE",
                                        start_date=start_date,
                                        end_date=end_date,
                                       )
    return npp_14c


def get_npp_flux


    store="radiocarbon_store",
                                species="dco2c14",
                                source="radd-eu",
                                domain="EUROPE",
                                start_date="2021-01-01",
                                end_date="2021-02-01",