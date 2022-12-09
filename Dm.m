function [D] = Dm(k)
D = zeros(6);
DAAs = 0;
DBBs = 0;
DBAs = 0;
DABs = 0;
for i = 1:3
    DAAs = DAAs+KAA(:,:,i)*exp(1j*(k* -SA(:,:,i))) + KAA(:,:,i+3)*exp(1j*(k* -SA(:,:,i+3)));
    DBBs = DBBs + KBB(:,:,i) * exp(1j * (k* -SB(:,:,i))) + KBB(:,:,i+3)*exp(1j*(k* -SB(:,:,i+3)));
    DABs = (DABs + KAB1(:,:,i) * exp(1j * (k* -FA(:,:,i))) + KAB3(:,:,i) * exp(1j * (k* -TA(:,:,i))) + ...
            KAB4(:,:,i) * exp(1j * (k* -LA(:,:,i))) + KAB4(:,:,i+3) * exp(1j * (k* -LA(:,:,i+3))));
    DBAs = (DBAs + KBA1(:,:,i) * exp(1j * (k* FA(:,:,i))) + KBA3(:,:,i) * exp(1j * (k* TA(:,:,i))) +  ...
            KBA4(:,:,i) * exp(1j * (k* LA(:,:,i))) + KBA4(:,:,i+3) * exp(1j * (k* LA(:,:,i+3))));
D(1:3, 4:6) = -DABs/sqrt(MN)/sqrt(MB);
D(4:6, 1:3) = -DBAs/sqrt(MN)/sqrt(MB);
D(1:3, 1:3) = (sum(KAB1,3)+sum(KAA,3)+sum(KAB3,3)+sum(KAB4f,3)+sum(KAB4s,3)-DAAs)/MB;
D(4:6, 4:6) = (sum(KBA1,3) + sum(KBB,3) + sum(KBA3,3) + sum(KBA4f,3) + sum(KBA4s,3) - DBBs)/MN;
end
end

