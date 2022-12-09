function [R,RB] = rotate(theta, r)
%         theta: ��ת�Ƕ�
%         r: λ���������K����
%         ��ת֮��ͨ����ת����U���Խ���A��ԭ��Ϊ���ĵĽ���ԭ�ӵ���Ӧ��K����ת������B��ԭ��Ϊ���ĵ�
%         theta: ��ת�Ƕ�
%         r: λ���������K����
%         ��ת֮��ͨ����ת����U���Խ���A��ԭ��Ϊ���ĵĽ���ԭ�ӵ���Ӧ��K����ת������B��ԭ��Ϊ���ĵ�
    R = []
    RB = []
    Um = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1]
    Umn = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1]
    U = [-1, 0, 0; 0, -1, 0; 0, 0, 1]
    for i = 1:fix(2*pi / theta)    % ȷ����ת����
        if size(r) == [3,1]              % ͨ���ж���ȷ����λ�����껹��K����
            r = Umn*r   % λ��������ת
            rb = r*(-1)                            % ����ԭ�Ӹı���λ��
        else
            r =Umn*r*Um             % ��K����Ĳ�����ͬλ������Ĳ�������
            rb =U*r*U 
        end
        R(:,:,i) =r                      % ��һ���б����棬������ֱ����ͣ����㷽��?
        RB(:,:,i)=rb                  %����Bԭ��Ϊ���ĵ�λ�ú�K����
    end
end



