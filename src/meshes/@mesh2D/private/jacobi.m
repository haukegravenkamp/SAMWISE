function jac = jacobi(a,b,n,x)

jac=ones(size(x));

if(n>0)
    
    al = 0.5*(a+b+2); be= (b-a)/(a+b+2);
    p  = al*(x-be);
    
    if(n>1)
        pm1=p;
        al = (2+a+b+1)*(2+a+b+2)/(4*(a+b+2) );
        be = (b^2-a^2)/((2+a+b)*(2+a+b+2));
        ga = (1+a)*(1+b)*(2+a+b+2)/(2*(a+b+2)*(2+a+b));
        
        p  = al*(x-be).*pm1 - ga;
        
        if(n>2)
            
            for i=2:n-1
                pm2=pm1;
                pm1=p;
                
                al = (2*i+a+b+1)*(2*i+a+b+2)/(2*(i+1)*(i+a+b+1) );
                be = (b^2-a^2)/((2*i+a+b)*(2*i+a+b+2));
                ga = (i+a)*(i+b)*(2*i+a+b+2)/((i+1)*(i+a+b+1)*(2*i+a+b));
                p  = al*(x-be).*pm1 - ga*pm2;
            end
            
        end
    end
    jac=p;
end


end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
