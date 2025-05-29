function LFI=LinearFISimple(act,map,lcRad,n)
    index=round(map*12/pi);
    index=index(lcRad+1:lcRad+n,lcRad+1:lcRad+n);
    dadt=devTuning(act,lcRad,n,map);
    dev=dadt(index+1);
    varMap=squeeze(var(act));
    LFI=(dev.*dev)./varMap;
end