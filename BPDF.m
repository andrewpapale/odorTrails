function NR=BPDF(X)
% X noise image
parxek=30;
paryek=30;
parx=1;
pary=1;
border=10;
NI=padarray(X,[parx+parxek pary+paryek],'symmetric');
            [m,n]=size(NI);            
    for i=parx+parxek+1:m-(parx+parxek)
    for j=pary+paryek+1:n-(pary+paryek) 
          t=0;
     if (NI(i,j)==0 | NI(i,j)==255)
         n_v=1;
       while n_v > 0       
            array_1=NI(i-(parx+t):(i+parx+t),j-(pary+t):j+pary+t);   
            array_2=array_1(array_1>0 & array_1<255); 
                 t1=0;
                 t2=0;            
            min_2=min(array_2());
            max_2=max(array_2());            
            if((min_2-border)~=0)
             t1=1;  
            end
             if((max_2+border)~=255)
              t2=1;
             end          
             if(t1==1 && t2==1)
               array_3= array_1(array_1>0 & array_1<255);         
             end            
             if(t1==1 && t2==0)
                array_3=array_1(array_1>0);         
             end               
             if(t1==0 && t2==1)
               array_3=array_1(array_1<255);         
             end
             if(t1==0 && t2==0)
                array_3=array_1;         
             end           
           [uv1,~,idx1] = unique(array_3);
            t3 = accumarray(idx1(:),1);
            max_value = max(t3);
           z = find(t3 ==  max_value);
           [a b]=size(z);
           if(a>0)
                n_v=0;
                t=0;
           if(a>1)
             c=ceil(a/2);
             NI(i,j)=uv1(z(c,1),1);   
           else
           NI(i,j)=uv1(z(a,1),1);    
           end
           else
               t=t+1;
              n_v=1;  
           end
       end
     end
    end
    end
      NR=NI(parx+parxek+1:m-(parx+parxek),pary+paryek+1:n-(pary+paryek));      
        