# import math
import lake as lk
import numpy as np
import copy
# reload(lk)
# dir(lk)


########################################################################
# get_args can get the interface informaiton
# of the f2py-converted Fortran interface
def get_args(doc_string):
    # split doc lines for processing
    docLines = np.array(
        doc_string.splitlines(),
        dtype=str)

    # get the information of input variables for SUEWS_driver
    posInput = np.where(
        np.logical_or(
            docLines == 'Parameters', docLines == 'Returns'))
    varInputLines = docLines[posInput[0][0] + 2:posInput[0][1] - 1]
    varInputInfo = np.array([[xx.rstrip() for xx in x.split(':')]
                             for x in varInputLines])
    dict_InputInfo = {xx[0]: xx[1] for xx in varInputInfo}
    dict_InOutInfo = {xx[0]: xx[1] for xx in varInputInfo if 'in/out' in xx[1]}

    # get the information of output variables for SUEWS_driver
    posOutput = np.where(docLines == 'Returns')
    varOutputLines = docLines[posOutput[0][0] + 2:]
    varOutputInfo = np.array([[xx.rstrip() for xx in x.split(':')]
                              for x in varOutputLines])
    dict_OutputInfo = {xx[0]: xx[1] for xx in varOutputInfo}

    # pack in/out results:
    dict_inout = {
        # 'input' and 'output' are dict's that store variable information:
        # 1. intent: e.g., input, in/output
        # 2. dimension: e.g., (366,7)
        'input': dict_InputInfo,
        'output': dict_OutputInfo,
        # 'var_input' and 'var_output' are tuples,
        # that keep the order of arguments as in the Fortran subroutine
        'var_input': tuple(varInputInfo[:, 0]),
        'var_inout': tuple(dict_InOutInfo.keys()),
        'var_output': tuple(varOutputInfo[:, 0])}

    return dict_inout


# get the inout info about `lake`
dict_info_lake = get_args(lk.lakemain.__doc__)
list_var_input = dict_info_lake['var_input']
list_var_inout = dict_info_lake['var_inout']
list_var_output = dict_info_lake['var_output']


# high-level wrapper: lake_cal_tstep
def lake_cal_tstep(dict_state_start, dict_met_forcing_tstep):
    # use single dict as input for suews_cal_main
    dict_input = dict_state_start.copy()
    dict_input.update(dict_met_forcing_tstep)
    dict_input = {k: dict_input[k] for k in list_var_input}

    # main calculation:
    res_suews_tstep = lk.lakemain(**dict_input)

    # update state variables
    dict_state_end = dict_state_start.copy()
    dict_state_end.update({var: copy.copy(dict_input[var])
                           for var in list_var_inout})

    # pack output
    dict_output = {k: v for k, v in zip(
        list_var_output, res_suews_tstep)}

    return dict_state_end, dict_output


# use a dict to hold parameters
dict_state = {}

# time step: dtime
dtime = np.array([60.])
dict_state.update(dtime=np.array([60.]))
# input: constants
# forc_hgt = np.array([10.])
# forc_hgt_q = np.array([10.])   # observational height of humidity [m]
# forc_hgt_t = np.array([10.])  # observational height of temperature [m]
# forc_hgt_u = np.array([10.])  # observational height of wind [m]

dict_state.update(forc_hgt=np.array([10.]),
                  forc_hgt_q=np.array([10.]),
                  forc_hgt_t=np.array([10.]),
                  forc_hgt_u=np.array([10.]))

# latitude (radians) latitude of Nuozhadu, 21.77 degN
# lat=np.deg2rad([21.7724])
dict_state.update(lat=np.deg2rad([21.7724]))

# Vertical Discretization
z_lake = np.array([[0.05, 0.15, 0.3, 0.5, 0.75, 1.15, 1.65, 2.2, 3, 4, 5.2,
                    6.7, 8.5, 11, 14, 17.5, 22.5, 29.5, 38.5, 49.5, 63.5,
                    81.5, 104.5, 134.5, 173.5]])
dict_state.update(z_lake=z_lake)

# dz_lake = np.array([[0.1,0.1,0.2,0.2,0.3,0.5,0.5,0.6,1,1,
#                      1.4,1.6,2,3,3,4,6,8,10,12,16,20,26,34,44]])
dz_lake = np.insert(np.diff(z_lake), 0, .1, 1)
lakedepth = np.array([195.5])  # column lake depth (m)
do_capsnow = np.array([0], dtype=np.int32)

dict_state.update(dz_lake=dz_lake, lakedepth=lakedepth, do_capsnow=do_capsnow)

# initialization
h2osno = np.array([0.])
snowdp = np.array([0.])
snl = np.array([0], dtype=np.int32)
# layer depth for snow & soil (m)
z = np.array([[0.00, 0.00, 0.00, 0.00, 0.00, 1.05, 1.15, 1.25, 1.35, 1.45,
               1.55, 1.65, 1.75, 1.85, 1.95]])
# layer thickness for soil or snow (m)
dz = np.array([[0.00, 0.00, 0.00, 0.00, 0.00, 0.1, 0.1, 0.1, 0.1, 0.1,
                0.1, 0.1, 0.1, 0.1, 0.1]])
# interface level below a "z" level (m)
zi = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.1, 1.2, 1.3, 1.4, 1.5,
                1.6, 1.7, 1.8, 1.9, 2.0]])

dict_state.update(h2osno=h2osno, snowdp=snowdp,
                  snl=snl, z=z, dz=dz, zi=zi)

# volumetric soil water (0<=h2osoi_vol<=watsat)[m3/m3]
h2osoi_vol = np.array(
    [[0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1,
      0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]])
# liquid water (kg/m2) NB: wrong variable name?
# h2osoi_vol = np.array(
#     [[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
#       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]])
h2osoi_liq = np.array(
    [[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]])
# ice lens (kg/m2)
h2osoi_ice = np.array(
    [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
t_grnd = np.array([294.])
# soil (or snow) temperature (Kelvin)
t_soisno = np.array([[273.15, 273.15, 273.15, 273.15, 273.15,
                      284.15, 284.15, 284.15, 284.15, 284.15,
                      284.15, 284.15, 284.15, 284.15, 284.15]])
# t_lake = np.array([[288.15, 288.15, 288.15, 288.15, 288.15,
# 288.15, 288.15, 288.15, 288.15, 288.15]])
t_lake = np.array([[294, 294, 294, 294, 294, 294, 293.8, 293.8,
                    293.8, 293.8, 293.8, 293.8, 293.8, 293.8, 293.8, 293.8,
                    293.8, 293.3, 292.1, 291.5, 290.7, 289.7, 287.1, 285.5,
                    284.5]])

dict_state.update(h2osoi_vol=h2osoi_vol,
                  h2osoi_liq=h2osoi_liq,
                  h2osoi_ice=h2osoi_ice,
                  t_grnd=t_grnd,
                  t_soisno=t_soisno,
                  t_lake=t_lake)

# top level eddy conductivity from previous timestep (W/m.K)
# eddy conductivity is smaller and varies between 0.01 and 100
savedtke1 = np.array([0.00001])
# mass fraction of lake layer that is frozen
lake_icefrac = np.array([[0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0]])
sand = np.array([100.]).repeat(10)
clay = np.array([0.]).repeat(10)
watsat = 0.489 - 0.00126 * sand
csol = (2.128 * sand + 2.385 * clay) / (sand + clay) * 1000000.  # J/(m3 K)
tkm = (8.80 * sand + 2.92 * clay) / (sand + clay)
tkmg = tkm ** (1. - watsat)
bd = (1. - watsat) * 2700.
tkdry = (0.135 * bd + 64.7) / (2700. - 0.947 * bd)
tksatu = tkmg * 0.57**watsat
watsat = watsat.reshape(1, watsat.shape[0])
csol = csol.reshape(1, csol.shape[0])
tkmg = tkmg.reshape(1, tkmg.shape[0])
tkdry = tkdry.reshape(1, tkdry.shape[0])
tksatu = tksatu.reshape(1, tksatu.shape[0])

dict_state.update(savedtke1=savedtke1,
                  lake_icefrac=lake_icefrac,
                  sand=sand,
                  clay=clay,
                  csol=csol,
                  tkm=tkm,
                  tkmg=tkmg,
                  bd=bd,
                  tkdry=tkdry,
                  tksatu=tksatu,
                  watsat=watsat)


# input and inoutflow
a = np.loadtxt('data/input.txt')
# a = a.reshape(a.shape[0], -1, 1)
array_forcing = np.reshape(np.loadtxt('data/input.txt'), (-1, 1, 13, 1))

# yeara = array_forcing[:, 0]
# montha = array_forcing[:, 1]
# daya = array_forcing[:, 2]
# houra = array_forcing[:, 3]
# forc_ta = array_forcing[:, 4]
# forc_pbota = array_forcing[:, 5]
# forc_psrfa = array_forcing[:, 6]
# forc_qa = array_forcing[:, 7]
# forc_ua = array_forcing[:, 8]
# forc_va = array_forcing[:, 9]
# forc_lwrada = array_forcing[:, 10]
# preca = array_forcing[:, 11]
# sabga = array_forcing[:, 12]

output8 = []
output20 = []

dict_state.update(x_lake_icefrac=dict_state['lake_icefrac'][0, 0])

dict_state_start = dict_state.copy()
# each hour for year 2015
for i in range(0, 8760, 1):
    # year = yeara[i]
    # month = montha[i]
    # day = daya[i]
    # hour = houra[i]
    # forc_t = forc_ta[i]
    # forc_pbot = forc_pbota[i]
    # forc_psrf = forc_psrfa[i]
    # forc_q = forc_qa[i]
    # forc_u = forc_ua[i]
    # forc_v = forc_va[i]
    # forc_lwrad = forc_lwrada[i]
    # prec = preca[i]
    # sabgc = sabga[i]
    print dict_info_lake['var_input']

    # use dict to hold te input variables
    list_var = ['year', 'month', 'day', 'hour',
                'forc_t', 'forc_pbot', 'forc_psrf', 'forc_q',
                'forc_u', 'forc_v', 'forc_lwrad', 'prec', 'sabgc']

    # i = 12
    dict_forcing_tstep = {var: val for var, val
                          in zip(list_var, array_forcing[i, 0])}

    # the following 1 is different for each minute: solar radiation * albedo
    dict_forcing_tstep.update(
        sabg=dict_forcing_tstep['sabgc'] * (
            1 - (0.6 * dict_state['x_lake_icefrac'])
            - ((1.0 - dict_state['x_lake_icefrac']) * 0.08)))

    dict_input_tstep = dict_forcing_tstep.copy()
    dict_input_tstep.update(dict_state)
    dict_input_tstep = {k: dict_input_tstep[k]
                        for k in list(dict_info_lake['var_input'])}
    # simulation for each minute
    for k in range(0, 60, 1):
        print i, k
        dict_state_start = {k: dict_state_start[k]
                            for k in list(dict_info_lake['var_inout'])}
        res_tstep = lk.lakemain(**dict_input_tstep)
        dict_state_end = dict_state_start.copy()
        dict_state_end.update({var: copy.copy(dict_input_tstep[var])
                               for var in dict_info_lake['var_inout']})


# TODO: use pandas to do the post-processing


    # kme, eflx_lwrad_net, eflx_gnet, eflx_sh_tot, eflx_lh_tot, t_ref2m, q_ref2m, taux, tauy, ram1, z0mg, z0hg, z0qg = lk.lakemain(
    #     watsat, tksatu, tkmg, tkdry, csol, forc_t, forc_pbot, forc_psrf, forc_hgt, forc_hgt_q, forc_hgt_t, forc_hgt_u, forc_q, forc_u, forc_v, forc_lwrad, prec, sabg, lat, z_lake, dz_lake, lakedepth, do_capsnow, h2osno, snowdp, snl, z, dz, zi, h2osoi_vol, h2osoi_liq, h2osoi_ice, t_grnd, t_soisno, t_lake, savedtke1, lake_icefrac)

#         # output at 8:00, which is the same as measured lake temperature
#         if hour == 7 and k == 59:
#             output8.append([year, month, day, eflx_gnet, eflx_lwrad_net, eflx_sh_tot, eflx_lh_tot, t_ref2m, q_ref2m, taux, tauy, ram1, z0mg, z0hg, z0qg,
#                             t_grnd[0], t_lake[0, 0], t_lake[0, 1], t_lake[0, 2], t_lake[0, 3], t_lake[0,
#                                                                                                       4], t_lake[0, 5], t_lake[0, 6], t_lake[0, 7], t_lake[0, 8], t_lake[0, 9],
#                             t_lake[0, 10], t_lake[0, 11], t_lake[0, 12], t_lake[0, 13], t_lake[0,
#                                                                                                14], t_lake[0, 15], t_lake[0, 16], t_lake[0, 17], t_lake[0, 18], t_lake[0, 19],
#                             t_lake[0, 20], t_lake[0, 21], t_lake[0,
#                                                                  22], t_lake[0, 23], t_lake[0, 24],
#                             snl, lake_icefrac[0, 0], snowdp, sabg,
#                             kme[0, 0], kme[0, 1], kme[0, 2], kme[0, 3], kme[0,
#                                                                             4], kme[0, 5], kme[0, 6], kme[0, 7], kme[0, 8], kme[0, 9],
#                             kme[0, 10], kme[0, 11], kme[0, 12], kme[0, 13], kme[0,
#                                                                                 14], kme[0, 15], kme[0, 16], kme[0, 17], kme[0, 18], kme[0, 19],
#                             kme[0, 20], kme[0, 21], kme[0, 22], kme[0, 23], kme[0, 24]])
#
#         # output at 20:00, which is the same as measured lake temperature
#         if hour == 19 and k == 59:
#             output20.append([year, month, day, eflx_gnet, eflx_lwrad_net, eflx_sh_tot, eflx_lh_tot, t_ref2m, q_ref2m, taux, tauy, ram1, z0mg, z0hg, z0qg,
#                              t_grnd[0], t_lake[0, 0], t_lake[0, 1], t_lake[0, 2], t_lake[0, 3], t_lake[0,
#                                                                                                        4], t_lake[0, 5], t_lake[0, 6], t_lake[0, 7], t_lake[0, 8], t_lake[0, 9],
#                              t_lake[0, 10], t_lake[0, 11], t_lake[0, 12], t_lake[0, 13], t_lake[0,
#                                                                                                 14], t_lake[0, 15], t_lake[0, 16], t_lake[0, 17], t_lake[0, 18], t_lake[0, 19],
#                              t_lake[0, 20], t_lake[0, 21], t_lake[0,
#                                                                   22], t_lake[0, 23], t_lake[0, 24],
#                              snl, lake_icefrac[0, 0], snowdp, sabg,
#                              kme[0, 0], kme[0, 1], kme[0, 2], kme[0, 3], kme[0,
#                                                                              4], kme[0, 5], kme[0, 6], kme[0, 7], kme[0, 8], kme[0, 9],
#                              kme[0, 10], kme[0, 11], kme[0, 12], kme[0, 13], kme[0,
#                                                                                  14], kme[0, 15], kme[0, 16], kme[0, 17], kme[0, 18], kme[0, 19],
#                              kme[0, 20], kme[0, 21], kme[0, 22], kme[0, 23], kme[0, 24]])
#
#
# # write output in txt named 'output8.txt'
# np.savetxt('output8.txt', output8, ["%4d", "%2d", "%2d", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.4f", "%10.4f", "%10.2f", "%10.6f", "%10.6f", "%10.6f",
#                                     "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f",
#                                     "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f",
#                                     "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f",
#                                     "%10.2f", "%10.2f", "%10.2f", "%10.2f",
#                                     "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f",
#                                     "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f",
#                                     "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f"], delimiter='  ', newline='\n')
#
# # write output in txt named 'output20.txt'
# np.savetxt('output20.txt', output20, ["%4d", "%2d", "%2d", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.4f", "%10.4f", "%10.2f", "%10.6f", "%10.6f", "%10.6f",
#                                       "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f",
#                                       "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f",
#                                       "%10.2f", "%10.2f", "%10.2f", "%10.2f", "%10.2f",
#                                       "%10.2f", "%10.2f", "%10.2f", "%10.2f",
#                                       "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f",
#                                       "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f",
#                                       "%10.8f", "%10.8f", "%10.8f", "%10.8f", "%10.8f"], delimiter='  ', newline='\n')
