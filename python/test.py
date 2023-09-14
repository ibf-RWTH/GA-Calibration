
import numpy as np

def get_material_varbound(material_ids: list) -> np.ndarray:
    """
    return varbound for calibration based on material ids
    """
    varbound = []
    for mat_id in material_ids:
        if mat_id in MatID2MatProps:
            mat_props_bounds = MatID2MatPropsBound[mat_id]
            for value in mat_props_bounds.values():
                varbound.append(value)
    varbound = np.array(varbound)
    return varbound# sample for all possible parameters
MatID2MatProps = {
        1: ['pw_fl','brgvec','v0','rho_0','km1','km2','c3'], 
        2: ['pw_fl','hdrt_0','crss_0','crss_s','pw_hd','Adir','Adyn','gam_c','pw_irr'],
        3: ['pw_fl','shrt_0', 'crss_0','k','crss_s','pw_hd','Adir','Adyn','gam_c','pw_irr'],
        4: ['pw_fl','shrt_0', 'crss_0','k','crss_s','pw_hd','Adir','Adyn','gam_c','pw_irr']
    }
# these are default bound    
MatProp2Bound = {
        'shrt_0': [0.001, 1000],
        'pw_fl': [1, 100],
        'hdrt_0': [800, 1200],
        'crss_0': [10, 190],
        'crss_s': [100, 500],
        'pw_hd' : [1, 3],
        'Adir': [0, 700],
        'Adyn': [0, 300]
}
# choose phases and the varrying parameters for each phase here
material_id = [2, 3]
MatID2MatProps[2] = ['pw_fl','hdrt_0','crss_0','crss_s']
MatID2MatProps[3] = ['pw_fl','hdrt_0','crss_0','crss_s','pw_hd' ]

MatID2MatPropsBound = {}
for mat_id in material_id:
    MatID2MatPropsBound[mat_id] = {}

for mat_id in material_id:
    if mat_id in MatID2MatProps:
        mat_props = MatID2MatProps[mat_id]
        for mat_prop in mat_props:
            MatID2MatPropsBound[mat_id][mat_prop] = MatProp2Bound[mat_prop]

MatID2MatPropsBound[2]['hdrt_0'] = [800, 1199]
varbound = get_material_varbound(material_id)
print(varbound)