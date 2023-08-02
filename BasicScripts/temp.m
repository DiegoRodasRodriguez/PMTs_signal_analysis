

q1=U(:,1);q2=U(:,2);qsc=U(:,3);t1=U(:,5);t2=U(:,6);
Dq=q1-q2;Dt=t1-t2;Sq=q1+q2;St=t1+t2;

I=find(qsc>1000);
[q1,q2,qsc,t1,t2,Dq,Dt,St,Sq]=cuts(I,q1,q2,qsc,t1,t2,Dq,Dt,St,Sq);