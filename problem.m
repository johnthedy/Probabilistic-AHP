function [nAHP,AHPraw]=problem(nprob)
if nprob==1
    nAHP = 4;
    AHPraw(1,:)={0 [5,7.1/9,1/8,1/7,3,3,4,7,8.7,8,5,6,7,8,9,1/5,1/4,2,3] 	[9,7,1/6,1/5,1/4,4,1/5,1/4,4,5,2,3,1/3,1/2,1,4,5,6,2,3,5,6]	[7,5,3,4,5,5,1/7,1/6,6,7,5,6,1,2,3,3,4,5,6,7,8,9]};
    AHPraw(2,:)={0 0 [6,4,4,5,3,1/7,1/6,1/3,1/2,1/6,1/5,1/8,1/7,1/7,1/6,1/5,7,8,2,3]	[1/7,1/3,8,9,4,1/8,1/7,1/5,1/4,1/5,1/4,1/5,1/4,1/3,1/6,1/5,1/4,8,9,5,6]};
    AHPraw(3,:)={0 0 0 [1/8,1/4,4,5,6,2,1/5,1/4,1,2,3,4,7,8,9,1,2,3,4,5,4,5]};
elseif nprob==2
    nAHP = 6;
    AHPraw(1,:)={0 [2,5,4,5,3,2,3,4,1,2,3,4,1/3,1/2,1.4,5,4,5] [4,7,5,6,7,7,7,8,6,7,4,4,5,6,3,4,5,7,8,6] [4,8,7,8,9,3,1/4,1/3,5,6,6,7,6,7,5,6,7,1,2,9] [4,7,6,7,8,5,5,6,4,5,6,3,4,5,4,5,6,2,3,7] [6,8,6,7,2,6,7,2,3,6,7,8,5,6,7,6,7,8]};
    AHPraw(2,:)={0 0 [5,8,6,7,8,3,6,7,7,8,7,8,5,6,7,8,9,6,7,5] [5,7,4,5,6,5,5,6,7,8,3,4,5,2,3,4,4,5,2] [5,6,5,6,7,3,3,4,6,7,5,3,4,4,5,6,1/6,1/5,3] [6,7,5,6,7,1/2,4,5,4,5,5,6,7,5,6,7,3,4,4]};
    AHPraw(3,:)={0 0 0 [1/3,1,1/4,1/3,3,1/3,1/2,2,3,1/5,1/4,1/3,1/2,1,1/7,1/6,1/5,1/4,1/3,1/4] [1/5,1/2,1/4,1/3,1/3,1/7,1/6,1/4,1/3,1/4,1/3,1/5,1/4,1/3,1/5,1/4,1/3,1/7,1/6,1/3] [1/3,1,1/6,1/5,1/4,1/5,1/4,1/5,1/4,1/3,1,2,3,1/3,1/2,1,1/5,1/4,1/2]};
    AHPraw(4,:)={0 0 0 0 [3,6,4,5,1/3,1/5,1/4,1/5,1/4,3,1/3,1/2,4,5,6,1/7,1/6,2] [2,5,2,3,1/5,1/4,1/3,1/6,1/5,5,1,2,3,5,6,7,1/2,1,3]};
    AHPraw(5,:)={0 0 0 0 0 [2,5,1/4,1/3,1/4,2,3,1/3,1/2,3,4,5,6,1,2,3,4,5,2]};
elseif nprob==3
    nAHP = 7;
    AHPraw(1,:)={0 [1/5,1/4,1/2,4,5,1/5,1/9,1/3,1/2,1/4,1/3]	[4,5,5,2,3,3,1/3,4,6,7]	[6,7,6,2,3,7,8,1/3,5,1/5]	[3,4,3,2,3,6,7,1/3,6,5]	[5,6,4,4,5,8,9,1/3,2,4]	[5,6,2,4,5,1/3,1/2,1/3,8,9,3]};
    AHPraw(2,:)={0 0 [7,6,1/6,1/5,5,6,9,5,6,5]	[7,7,1/6,1/5,7,8,9,6,7,1/3]	[7,4,1/4,1/3,6,7,9,8,3]	[7,5,1/4,1/3,9,9,3,4,4]	[7,3,1/4,1/3,2,3,9,9,4,5]};
    AHPraw(3,:)={0 0 0 [2,3,2,2,3,5,1/3,2,3,1/8]	[1/4,1/3,1/3,2,3,2,1/3,4,5,1/5]	[1/3,1/2,1/2,2,3,5,6,1/3,1/3,1/2,1/5,1/4]	[1/4,1/4,2,3,1/6,1/5,1/3,7,8,1/2]};
    AHPraw(4,:)={0 0 0 0 [1/6,1/5,1/4,2,3,1/2,1/3,3,4,7]	[1/4,1/3,3,4,5,2,3,1/3,1/5,1/4,5,6]	[1/5,1/4,1/5,2,3,4,1/7,1/3,6,7,5]};
    AHPraw(5,:)={0 0 0 0 0 [3,2,2,3,3,4,1/2,1/6,1/5,2]	[2,3,1/2,2,3,1/8,1/7,1/3,7,8,2,3]};
    AHPraw(6,:)={0 0 0 0 0 0 [1/3,1/3,1/4,1/3,1/9,1/8,1,7,8,2]};
end
end