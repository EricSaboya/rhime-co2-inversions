def rhime_co2_dictionary_keys():
    """
    -------------------------------------------------------
    Returns RHIME input dictionaries with descriptions
    of each variable. Dictionary keys are the input 
    variable and the corresponding entry is a description 
    of the input. 
    -------------------------------------------------------
    Returns
        - obs_dict_keys (dict)
        - flux_dict_keys (dict)
        - fp_dict_keys (dict)
        - bc_dict_keys (dict)
        - basis_dict_keys (dict)
        - mcmc_inputs_dict_keys (dict)
        - c14_dict_keys (dict)

    -------------------------------------------------------
    """
    # Observation dictionary keys
    obs_dict_keys = {"species":                            "Atmospheric trace gas species of interest (str)",
                     "site":                               "Measurement sites (list of str)",
                     "inlet":                              "Measurement site inlets (list of str)",
                     "averaging_period":                   "Measurement data averaging period (list of str)",
                     "instrument":                         "Instrument used for measurements (list of str)",
                     "data_level":                         "Data quality level (list of str) - typically used only for ICOS data",
                     "store":                              "OpenGHG object store name (str)",
                     "calibration_scale":                  "Measurement data calibration scale (str)",
                     "start_date":                         "Date for the start of period of interest (str) e.g. '2022-01-01' ",
                     "end_date":                           "Date for the end of period of interest (str) e.g. '2022-02-01' ",
                     "filters":                            "Measurement data filters to apply",
                    }
    
    # Flux dictionary keys
    flux_dict_keys = {"species":                           "Atmospheric trace gas species of interest (str)",
                      "domain":                            "Model domain name (str)",
                      "source":                            "Flux source name (list of str)",
                      "start_date":                        "Date for the start of period of interest (str) e.g. '2022-01-01' ",
                      "end_date":                          "Date for the end of period of interest (str) e.g. '2022-02-01' ",
                      "store":                             "OpenGHG object store name (str)",
                      "flux_sf":                           "Dictionary (keys = source) of multiplicative scale factors to apply to each flux sector",
                      "emisource_to_sector":               "Dictionary aligning flux sources with sector names e.g. {'fossil': 'edgar-fossil'}",
                     }

    # Footprint dictionary keys
    fp_dict_keys = {"species":                             "Atmospheric trace gas species of interest (str)",
                    "domain":                              "Model domain name (str)",
                    "site":                                "Measurement sites (list of str)",
                    "fp_height":                           "Footprint model inlet height (list of str)",
                    "start_date":                          "Date for the start of period of interest (str) e.g. '2022-01-01' ",
                    "end_date":                            "Date for the end of period of interest (str) e.g. '2022-02-01' ",
                    "store":                               "OpenGHG object store name (str)",
                   }
    
    # Boundary conditions dictionary keys
    bc_dict_keys = {"species":                             "Atmospheric trace gas species of interest (str)",
                    "domain":                              "Model domain name (str)",
                    "bc_input":                            "Source name of the BC inputs",
                    "bc_freq":                             "Frequency over which to solve the BC values",
                    "start_date":                          "Date for the start of period of interest (str) e.g. '2022-01-01' ",
                    "end_date":                            "Date for the end of period of interest (str) e.g. '2022-02-01' ",
                    "store":                               "OpenGHG object store name (str)",
                    "bc_sf":                               "A priori multiplicative scaling factor to apply to BCs. Can be str or dict. with value for each boundary",
                   }

    # Basis functions dictionary keys
    basis_dict_keys = {"fp_basis_case":                    "Name of basis function file to use (str, optional)",
                       "basis_directory":                  "Directory of where basis function file is kept (str, optional)",
                       "fp_basis_algorithm":               "Name of algorithm to apply to flux*fps (str)",
                       "nbasis":                           "No. of basis functions per sector e.g. [40, 30, 10] or 40",
                       "bc_basis_case":                    "Name of basis function file to use for BCs",
                       "bc_basis_directory":               "Name of basis function directory",
                      }

    # MCMC dict keys
    mcmc_inputs_dict_keys = {"xprior":                     "Dictionary of dictionaries with prior PDF params for each flux sector",
                             "bcprior":                    "Dictionary of prior PDF params for BC flux sector",
                             "sigprior":                   "Dictionary of model error hyperparameter scaling PDF",
                             "add_offset":                 "Option to solve for an additive bias in the inversions (bool)",
                             "offsetprior":                "Dictionary of offset prior PDF params (if needed)",
                             "nit":                        "No. of MCMC iterations (int)",
                             "burn":                       "No. of MCMC iterations to discard (int)",
                             "tune":                       "No. of MCMC iterations for tuning (int)",
                             "nchain":                     "No. of MCMC chains to run (int)",
                             "sigma_per_site":             "Option to solve model error for each site (bool)",
                            }
                        
    # Radiocarbon dictionary keys
    c14_dict_keys = {"c14_resp_sig": {"additive_bias":     "Additive bias to apply to the respiration c14 signature",
                                     },
            
                     "c14_gpp_sig": {"additive_bias":      "Additive bias to apply to the GPP c14 signature",
                                    },
            
                     "c14_bg": {"additive_bias":           "Additive bias to apply to the basckground c14 signature",,
                               },
                     
                     "c14_nuclear" : {"source": "", 
                                     }
                    }


def paths_to_check():
    print("---- calculate_basis_functions.py ----")
    print("load_landsea_indices() --> path to UKMO land-sea definitions grid")
    print("bucketbasisfunction() --> default path to outputdir")
    print("basis() --> openghginv_path")
    print("basis_boundary_conditions() --> openghginv_path")
    print("")
    print("---- data_checks.py ----")
    print("Make sure 'pint' is installed!")
    print("")
    print("---- example_run_rhime.py ----")
    print("sys.path.append --> change to local dir of rhime-co2-inversions")
    print("main() --> outputpath where inversion outputs are saved")
    print("main() --> coutryfile name where PARIS country definitions file is saved")
    print("")
    print("---- get_co2_data.py ----")
    print("sys.path.append --> change to local dir of rhime-co2-inversions")
    print("")
    print("---- rhime_co2.py ----")
    print("sys.path.append --> change to local dir of rhime-co2-inversions")
    print("")
    print("---- sensitivity.py ----")
    print("sys.path.append --> change to local dir of rhime-co2-inversions")
    print("")
    
    