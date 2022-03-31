function [Acc,Sen,Spe]=ConMax(labelTargetTest,label)
aa=(label-labelTargetTest);aa=aa(aa==0); Acc=(length(aa)/length(label))*100;
 CMax=confusionmat(labelTargetTest,label);
 if length(CMax)==2
     Tp=CMax(1,1);Tn=CMax(2,2);Fn=CMax(2,1);Fp=CMax(1,2);
     Sen=(Tp/(Tp+Fn))*100;Spe=(Tn/(Tn+Fp)*100);
 end
end