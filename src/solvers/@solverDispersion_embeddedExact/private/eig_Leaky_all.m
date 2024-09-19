function [lambda, tmp_lambda, X] = eig_Leaky_all(L2,L1,L0,M,R,typeCoupling,kappa,whn,dofD)

switch typeCoupling

    case 'F' % fluid on one side or the same fluid on both sides
        R = R{1};
        kfh1 = kappa;
        if nargout>2
            if any(dofD,'all')
                [lambda, tmp_lambda,X] = eig_Leaky_single_linear(L2,L0,M,R,kfh1,whn,dofD);
            else
                [lambda, tmp_lambda,X] = eig_Leaky_single(L2,L1,L0,M,R,kfh1,whn);
            end
        else
            if any(dofD,'all')
                [lambda, tmp_lambda] = eig_Leaky_single_linear(L2,L0,M,R,kfh1,whn,dofD);
            else
                [lambda, tmp_lambda] = eig_Leaky_single(L2,L1,L0,M,R,kfh1,whn);
            end
        end

    case 'FF'  % two fluids
        R1 = R{1};
        R2 = R{2};
        kfh1 = kappa(1);
        kfh2 = kappa(2);
        if nargout>2
            [lambda, tmp_lambda, X] = eig_Leaky(L2,L1,L0,M,R1,R2,kfh1,kfh2,whn);
        else
            [lambda, tmp_lambda] = eig_Leaky(L2,L1,L0,M,R1,R2,kfh1,kfh2,whn);
        end

    case 'S'  % solid on one side
        R1 = R{1};
        R2 = R{2};
        kappal = kappa(1);
        kappat = kappa(2);
        opts.refine = 0;
        if nargout>2
            [lambda, tmp_lambda, X] = eig_LeakySolid(L2,L1,L0,M,R1,R2,kappal,kappat,whn,opts);
        else
            [lambda, tmp_lambda] = eig_LeakySolid(L2,L1,L0,M,R1,R2,kappal,kappat,whn,opts);
        end

    case 'SS' % solid on both sides
        kappaL1 = kappa(1);
        kappaT1 = kappa(2);
        kappaL2 = kappa(3);
        kappaT2 = kappa(4);
        opts.refine = 0;
        opts.rand_orth = 1;
        opts.showgap = 0;
        opts.solver = 'eig';
        opts.twosideRQ = 1;
        if nargout>2
            [lambda, tmp_lambda, X] = eig_LeakySolidSolid(L2,L1,L0,M,R,kappaL1,kappaT1,kappaL2,kappaT2,whn,opts);
        else
            [lambda, tmp_lambda] = eig_LeakySolidSolid(L2,L1,L0,M,R,kappaL1,kappaT1,kappaL2,kappaT2,whn,opts);
        end

    case 'FS' % one side solid, one side fluid
        Rp = R{1};
        Rs = R{2};
        Rf = R{3};
        kappaL1 = kappa(1);
        kappaT1 = kappa(2);
        kappaF = kappa(3);
        opts.rand_orth = 1;
        opts.refine = 0;
        opts.showgap = 0;
        opts.twosideRQ = 1;
        opts.maxlogres = 1e-2;
        if nargout>2
            [lambda, tmp_lambda, X] = eig_LeakySolidFluid(L2,L1,L0,M,Rp,Rs,Rf,kappaL1,kappaT1,kappaF,whn,opts);
        else
            [lambda, tmp_lambda   ] = eig_LeakySolidFluid(L2,L1,L0,M,Rp,Rs,Rf,kappaL1,kappaT1,kappaF,whn,opts);
        end

end

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
