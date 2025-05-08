output_path = 'C:\Users\smirhe1\OneDrive - Louisiana State University\ce7701\1D_model_flow';
name_output = '360x1';
trim = vs_use([output_path '\trim-' name_output '.dat'],'quiet');
bed_level=vs_let(trim,'map-sed-series','DPS','quiet');
Z0 = bed_level(1,2,:);
Z0 = squeeze(bed_level(1,2,:));
plot(Z0)

