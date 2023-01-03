function [hipFlex, kneeFlex] = getJointAngles(hip,knee,ankle)
v1 = knee - hip;
v2 = ankle - knee;

hipFlex = atand(v1(1)/v1(2));
% kneeFlex = atand(v2(1)/v2(3));
kneeFlex = acosd(v2*v1'/(norm(v2)*norm(v1)));


end