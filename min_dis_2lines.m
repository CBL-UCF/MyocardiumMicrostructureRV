function dis_min = min_dis_2lines(l1p1,l1p2,l2p1,l2p2)

%To find the minimum distance between two lines connecting 4
%points

%Variable name pattern: l1p2 = line 1 point 2

%Consider line 1 and line 2 in in parametric form l1p1+s1v1 and l2p1+s2v2
%(v1 and v2 are vectors)

%v= (v1 crossproducr v2) is dirrection vector for line M which is
%orthogonal to both line1 and line2

%get projections of points l1p1(m1) and l2p1(m2) into M line and calculate the
%distance between those two points will give the minimum distance between
%two lines



v1=l1p2-l1p1;
v2=l2p2-l2p1;

mag_v1=norm(v1);
mag_v2=norm(v2);

U_v1 = v1/mag_v1;
U_v2 = v2/mag_v2;

vcross = cross(v1,v2);

%make closest point to the origin is first point of the line
if norm(l1p1)<norm(l1p2)
    v1_1=l1p1;
    v1_2=l1p2;
else
    v1_1=l1p2;
    v1_2=l1p1;
end

if norm(l2p1)<norm(l2p2)
    v2_1=l2p1;
    v2_2=l2p2;
else
    v2_1=l2p2;
    v2_2=l2p1;
end
%If two lines are parallel,check weather they overlap
v1 = v1_2-v1_1;
v2 = v2_2-v2_1;

denom=norm(vcross)^2;

if denom==0   %Two lines are parallel condition
    d1=dot(U_v1,v1_1-v2_1);
    d2=dot(U_v1,v1_1-v2_2);
    
    if d1<=0 && d2<=0  %line 2 before line 1 condition
        if abs(d1)<abs(d2)
            dis_min = norm(v1_1-v2_1);
        else
            dis_min = norm(v1_1-v2_2);
        end 
        
    elseif d1>=mag_v1 && d2>=mag_v2  %line2 after line1 condition
        if abs(d1)<abs(d2)
            dis_min = norm(v1_2-v2_1);
        else
            dis_min = norm(v1_2-v2_2);
        end
    else  %two line segments overlap condition
        v3=v1_1+d1*U_v1-v2_1;
        dis_min = norm(v3);
    end
    
else %Two lines are not parallel condition
    %find unit vectors to express any point in v1 in the form of v1_1+t1*Unit(V1)
    
    uu =  v2_1-v1_1;
    
    det1 = det([uu;U_v2;vcross]);
    det2 = det([uu;U_v1;vcross]);
    
    % Used reference - Sheetlet model/Resources
    v22 = dot(v2,v2);
    v11 = dot(v1,v1);
    v21 = dot(v2,v1);
    v21_1 = dot(uu,v1);
    v21_2 = dot(uu,v2);
    denom1 = v21*v21 - v22*v11;
    
    if isalmost(denom,0.0, 1e-4)
        t1=0;
        t2=(v11*t1 - v21_1)/v21;
    else
        t1 = (v21_2*v21-v22*v21_1)/denom1;
        t2 = (-v21_1*v21+v11*v21_2)/denom1;
    end
    
    t1 = max(min(t1,1),0);
    t2 = max(min(t2,1),0);
    
    t1=t1*mag_v1;
    t2=t2*mag_v2;
        


    
    %find closest points on two lines m1,m2
    %t1 = dot(cross(U_v2,vcross),(-v1_1+v2_1))/dot(vcross,vcross); %References in model folder
    %t1 = det1/denom
    m1 = v1_1+t1*U_v1;
    
    %t2 = dot(cross(U_v1,vcross),(-v1_1+v2_1))/dot(vcross,vcross);  %References in model folder
    %t2 = mag_v2*det2/denom
    m2 = v2_1+t2*U_v2;
    
    %Check whether m1 and m2 on line segments
        if dot((v1_1-m1),(v1_2-m1))<0 && dot((v2_1-m2),(v2_2-m2))<0
            %both points on line segments condition
            dis_min = norm(m1-m2);
        elseif dot((v1_1-m1),(v1_2-m1))>0 && dot((v2_1-m2),(v2_2-m2))<0
            %Only m2 on linesegment condition
            if dot((v1_1-m2),U_v1)<dot((v1_2-m2),U_v1)
                dis_min = norm(v1_1-m2);
            else
                dis_min = norm(v1_2-m2);
            end
        elseif dot((v1_1-m1),(v1_2-m1))<0 && dot((v2_1-m2),(v2_2-m2))>0
            %Only m1 on line segement condition
            if dot((v2_1-m1),U_v2)<dot((v2_2-m1),U_v2)
                dis_min = norm(v2_1-m1);
            else
                dis_min = norm(v2_2-m1);
            end
        else
            %Both points not on line segment condition
            if t1<0
                m1=v1_1;
            elseif t1>mag_v1
                m1=v1_2;
            end
            if t2<0
                m2=v2_1;
            elseif t2>mag_v2
                m2=v2_2;
            end
            dis_min = norm(m1-m2);
        end
    end
end



