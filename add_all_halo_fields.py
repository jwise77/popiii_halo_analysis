import yt
import ytree
import numpy as np
from yt.utilities.physical_constants import planck_constant_cgs
import halo_analysis as ha
from file_locations import *
yt.enable_parallelism()

# Add J21_LW and J_Lyman fields to yt
sigma_H2I = yt.YTQuantity(3.71e-18, 'cm**2')
sigma_HI = yt.YTQuantity(6.3e-18, 'cm**2')
J21_norm = yt.YTQuantity(1e21, 'cm**2/erg')
E_LW = yt.YTQuantity(12.4, 'eV')
E_HI = yt.YTQuantity(13.6, 'eV')
nu_H = yt.YTQuantity(3.2881e15, 'Hz')
default_central_radius = yt.YTQuantity(5, 'pc')

def _J21_LW(field, data):
    return J21_norm*data['enzo', 'H2I_kdiss']*E_LW/(4.0*np.pi*np.pi*sigma_H2I \
                                                    *nu_H)
yt.add_field(name=('gas', 'J21_LW'), function=_J21_LW, units='dimensionless', \
             display_name='J$_{21}$', sampling_type='cell')

def _J_Lyman(field, data):
    return data['enzo', 'HI_kph']*planck_constant_cgs/(4.0*np.pi*np.pi*sigma_HI)

yt.add_field(name=('gas', 'J_Lyman'), function=_J_Lyman, units='erg/cm**2', \
             display_name='J(Lyman)', sampling_type='cell')

# Load ytree arbor to update.
a = ytree.load(tree_file)
# Load full timeseries of datasets.
ts = yt.load(param_files)

# Add all analysis fields to the arbor.

a.add_analysis_field('cosmological_time', units='yr')

# Full sphere fields
# Mass Fields
a.add_analysis_field('sphere_mass', units='Msun')
a.add_analysis_field('sphere_dark_matter_mass', units='Msun')
a.add_analysis_field('sphere_gas_mass', units='Msun')
a.add_analysis_field('sphere_stellar_mass', units='Msun')
a.add_analysis_field('sphere_popiii_mass', units='Msun')
a.add_analysis_field('mass_threshold', units='dimensionless')
# Growth history fields
a.add_analysis_field('sphere_mass_growth', units='1/yr')
a.add_analysis_field('sphere_mass_growth_derivative', units='1/yr**2')
a.add_analysis_field('number_of_mergers', units='dimensionless')
# Pop III particle fields
a.add_analysis_field('sphere_popiii_number', units='dimensionless')
a.add_analysis_field('sphere_popiii_remnants', units='dimensionless')
# Metallicity fields
a.add_analysis_field('sphere_gas_metal_fraction', units='dimensionless')
a.add_analysis_field('sphere_stellar_metal_fraction', units='dimensionless')
# Gas inflow/outflow fields
a.add_analysis_field('sphere_inflow_mass', units='Msun/yr')
a.add_analysis_field('sphere_outflow_mass', units='Msun/yr')
a.add_analysis_field('sphere_inflow_metal_fraction', units='dimensionless')
a.add_analysis_field('sphere_outflow_metal_fraction', units='dimensionless')
# Radiation fields
a.add_analysis_field('sphere_J21_LW', units='dimensionless/(s*Hz)')
a.add_analysis_field('sphere_J_Lyman', units='erg/(cm**2*s*Hz)')

# Half radius sphere fields
#Mass fields
a.add_analysis_field('half_sphere_mass', units='Msun')
a.add_analysis_field('half_sphere_dark_matter_mass', units='Msun')
a.add_analysis_field('half_sphere_gas_mass', units='Msun')
a.add_analysis_field('half_sphere_stellar_mass', units='Msun')
a.add_analysis_field('half_sphere_popiii_mass', units='Msun')
# Pop III particle fields
a.add_analysis_field('half_sphere_popiii_number', units='dimensionless')
a.add_analysis_field('half_sphere_popiii_remnants', units='dimensionless')
# Metallicity fields
a.add_analysis_field('half_sphere_gas_metal_fraction', \
                     units='dimensionless')
a.add_analysis_field('half_sphere_stellar_metal_fraction', \
                     units='dimensionless')
# Gas inflow/outflow fields
a.add_analysis_field('half_sphere_inflow_mass', units='Msun/yr')
a.add_analysis_field('half_sphere_outflow_mass', units='Msun/yr')
a.add_analysis_field('half_sphere_inflow_metal_fraction', \
                     units='dimensionless')
a.add_analysis_field('half_sphere_outflow_metal_fraction', \
                     units='dimensionless')
# Radiation fields
a.add_analysis_field('half_sphere_J21_LW', units='dimensionless/(s*Hz)')
a.add_analysis_field('half_sphere_J_Lyman', units='erg/(cm**2*s*Hz)')

# Central propety fields
a.add_analysis_field('central_density', units='g/cm**3')
a.add_analysis_field('central_H2_fraction', units='dimensionless')

# Statistics
a.add_analysis_field('number_of_cells', units='dimensionless')
a.add_analysis_field('truncated_center_data', units='dimensionless')


# Load array of redshifts for each timestep.
ytree_redshifts = np.load('%s/ytree_redshifts.npy' % (storage_dir))

# Loop through each redshift and add all fields for all halos at that redshift.
start = len(ts) - ytree_redshifts.size
for i, redshift in enumerate(ytree_redshifts):

    # Set dataset corresponding to this redshift.
    ds = ts[start+i]
    # Select halos at this redshift
    hlist = a.select_halos("tree['tree', 'redshift'] == {}".format(redshift), \
                           fields=['redshift'])

    for halo in hlist:

        pos, rad = ha.halo_coords(ds, halo)

        sph = ds.sphere(pos, rad)
        halo['number_of_cells'] = sph['ones'].size
        # If the inner sphere (radius = rad/2) is resolved by less than 2 cells, flag it and increase the radius
        mindx = sph['dx'].min()
        if rad < 4*mindx:
            half_rad = 2*mindx
            halo['truncated_center_data'] = 1
        else:
            half_rad = rad/2
            halo['truncated_center_data'] = 0
        half_sph = ds.sphere(pos, half_rad)
        
        # Cosmological time is used for mass growth rate calculations so should be added early.
        halo['cosmological_time'] = ds.current_time.to('yr')
        
        # Central property fields
        inner_rad = max(default_central_radius, sph['dx'].min())
        halo['central_density'], halo['central_H2_fraction'] = ha.central_density_and_H2(ds, halo, inner_rad=inner_rad)

        # Full sphere fields
        # Mass, metallicity, inflow/outflow, and radiation fields
        ha.masses_and_metallicities(ds, sph, halo=halo)
        ha.inflow_outflow(ds, sph, pos, rad, halo=halo)
        ha.radiation(ds, sph, halo=halo)
        # Mass and radiation fields must be added before mass threshold field.
        halo['mass_threshold'] = ha.mass_threshold(halo)
        # Growth history fields
        halo['sphere_mass_growth'] = ha.specific_growth_rate(halo)
        halo['sphere_mass_growth_derivative'] = ha.mass_growth_derivative(halo)
        # The number of mergers is determined by the number of ancestor halos at the previous data output that have merged into the current halo.
        # If a halo has no ancestors, it is new and so had no recent mergers.
        try:
            halo['number_of_mergers'] = len(list(halo.ancestors)) - 1
        except:
            halo['number_of_mergers'] = 0
        
        # Half sphere fields
        # Mass, metallicity, inflow/outflow, and radiation fields.
        ha.masses_and_metallicities(ds, sph, halo=halo, is_half=True)
        ha.inflow_outflow(ds, sph, pos, half_rad, halo=halo, is_half=True)
        ha.radiation(ds, sph, halo=halo, is_half=True)
        # for f in a.field_list:
        #     print(f, halo[f])
        # import pdb; pdb.set_trace()
        # break

    # Save arbor with updated field values.
    # JHW note (12 Feb 2021): The looping needs to be changed, so that the tree objects can be given to save_arbor to save the analysis fields. See https://ytree.readthedocs.io/en/latest/Fields.html#id3
    a.save_arbor(trees=a[:], filename=full_arbor_file)
    # break
