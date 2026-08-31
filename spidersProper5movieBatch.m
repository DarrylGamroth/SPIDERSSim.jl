

waves = (800:250:1800)*1e-9;
for k=1:numel(waves);
    tmp = spidersProper5movie(waves(k), 1024 );
end
