"""M2 exploratory cerebral waterfall and coronary compression closures.

These are model-form sensitivities, not established minimum perfusion laws.
All flows still use the regional incidence matrix: a closed pressure window
stops flow, never deletes a tissue, reservoir or blood volume.
"""
import math

from .regional import RegionalModel


def waterfall(pa, pv, pc, resistance):
    """Signed flow [mL/s] for endpoint/closing pressures [mmHg].

    Below both endpoint pressures the threshold has no effect. Between them
    it limits the driving pressure. Above both endpoints it stops flow. The
    symmetric reverse branch preserves passivity and donor reversal.
    """
    if not all(math.isfinite(v) for v in (pa,pv,pc,resistance)) or resistance<=0:
        raise ValueError('Finite pressures and positive resistance required')
    return (max(pa-max(pv,pc),0)-max(pv-max(pa,pc),0))/resistance


class PressureRegionalModel(RegionalModel):
    """Same 15 states as RegionalModel; only three branch flow laws change."""
    def __init__(self, source_path, project_path, overlay_path, closure=None):
        super().__init__(source_path,project_path,overlay_path)
        brain=self.overlay['brain']; coronary=self.overlay['coronary']
        settings=dict(brain_enabled=True,coronary_enabled=True,
                      ICP_mmHg=brain['ICP_mmHg'],tone_offset_mmHg=brain['tone_offset_mmHg'],
                      brain_resistance_multiplier=1.,kappa=coronary['compression_kappa'])
        if closure:
            if set(closure)-set(settings): raise ValueError('Unknown closure option')
            settings.update(closure)
        for k in ('brain_enabled','coronary_enabled'):
            if not isinstance(settings[k],bool): raise ValueError('Closure enable flags must be boolean')
        for k in ('ICP_mmHg','tone_offset_mmHg','brain_resistance_multiplier','kappa'):
            if not math.isfinite(settings[k]) or settings[k]<0: raise ValueError('Invalid closure parameter')
        if settings['brain_resistance_multiplier']==0: raise ValueError('Brain resistance must be positive')
        anchor=next(r['target_flow_mL_s'] for r in self.project['branches'] if r['name']=='brain')
        # Calibrate ONCE at the nominal M2 closing pressure, using the revised
        # regional allocation. Sweeping ICP/tone must not recalibrate R and
        # thereby conceal the pressure sensitivity being investigated.
        nominal_pc=brain['ICP_mmHg']+brain['tone_offset_mmHg']
        nominal_r=(brain['resistance_anchor_AO_mmHg']-max(brain['resistance_anchor_VC_mmHg'],nominal_pc))/anchor
        self.brain_resistance=nominal_r*settings['brain_resistance_multiplier']
        if self.brain_resistance<=0: raise ValueError('Invalid cerebral calibration')
        self.closure_settings={**settings,'brain_nominal_resistance_mmHg_s_mL':nominal_r,
                               'brain_effective_resistance_mmHg_s_mL':self.brain_resistance}

    def pressures(self,t,state,active=True):
        p,m=super().pressures(t,state,active=active)
        p['PERI']=m['pericardial_pressure_mmHg']
        return p,m

    def closing_pressures(self,p):
        s=self.closure_settings
        return {'brain':s['ICP_mmHg']+s['tone_offset_mmHg'],
                **{name:max(p['PERI'],0)+s['kappa']*max(p[node]-p['PERI'],0)
                   for name,node in [('coronary_left','LV'),('coronary_right','RV')]}}

    def flows(self,t,state,open_valves,active=True,pressures=None):
        p=self.pressures(t,state,active=active)[0] if pressures is None else pressures
        q=super().flows(t,state,open_valves,active=active,pressures=p)
        pc=self.closing_pressures(p)
        if self.closure_settings['brain_enabled']:
            q['brain']=waterfall(p['AO'],p['VC'],pc['brain'],self.brain_resistance)
        if self.closure_settings['coronary_enabled']:
            for name in ('coronary_left','coronary_right'):
                q[name]=waterfall(p['AO'],p['VC'],pc[name],self.branch_by_name[name].resistance)
        return q
