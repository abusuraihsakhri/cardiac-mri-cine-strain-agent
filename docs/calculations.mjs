export function computeMetrics(input) {
  const n = (key) => {
    const value = Number(input[key]);
    if (!Number.isFinite(value)) throw new Error(`${key} must be a number.`);
    return value;
  };
  const edv = n('edv');
  const esv = n('esv');
  const hr = n('hr');
  const bsa = n('bsa');
  const mass = n('mass');
  if (edv <= 0 || esv <= 0 || edv <= esv) throw new Error('EDV must be greater than ESV, and both must be positive.');
  if (hr <= 0 || bsa <= 0 || mass < 0) throw new Error('Heart rate and BSA must be positive; LV mass cannot be negative.');
  const sv = edv - esv;
  const co = sv * hr / 1000;
  const result = {
    lvef: sv / edv * 100,
    sv,
    edvi: edv / bsa,
    esvi: esv / bsa,
    svi: sv / bsa,
    co,
    ci: co / bsa,
    lvmi: mass / bsa,
    gls: n('gls'),
    gcs: n('gcs'),
    grs: n('grs'),
    ecv: null,
  };
  const ecvKeys = ['nativeT1', 'postMyoT1', 'preBloodT1', 'postBloodT1', 'hct'];
  const supplied = ecvKeys.map((key) => String(input[key] ?? '').trim() !== '');
  if (supplied.some(Boolean) && !supplied.every(Boolean)) throw new Error('ECV requires all four T1 values and hematocrit.');
  if (supplied.every(Boolean)) {
    const nativeT1 = n('nativeT1');
    const postMyoT1 = n('postMyoT1');
    const preBloodT1 = n('preBloodT1');
    const postBloodT1 = n('postBloodT1');
    const hct = n('hct');
    if ([nativeT1, postMyoT1, preBloodT1, postBloodT1].some((value) => value <= 0)) throw new Error('T1 values must be positive.');
    if (hct <= 0 || hct >= 100) throw new Error('Hematocrit must be between 0 and 100%.');
    const deltaMyo = 1000 / postMyoT1 - 1000 / nativeT1;
    const deltaBlood = 1000 / postBloodT1 - 1000 / preBloodT1;
    if (deltaBlood <= 0) throw new Error('Post-contrast blood R1 must exceed pre-contrast blood R1.');
    result.ecv = (1 - hct / 100) * (deltaMyo / deltaBlood) * 100;
    if (result.ecv < 0) throw new Error('Calculated ECV is negative; check the paired T1 values.');
  }
  return result;
}
