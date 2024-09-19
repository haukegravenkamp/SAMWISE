function glb = fluidStructureCoupling(obj,glb,pdes,mat,intTab,opt)
% include coupling terms between elastic and acoustic layers in the FE-
% matrices. Only acts explicitly on the PDE classes pdeElasticity and
% pdeAcoustics. If other PDEs with their own interactions are included
% later, it should be easy to include them here.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
return

% TO DO

% 
% for i = 1:numel(glb)
%     glb(i) = getCouplingTerms(obj,glb(i),pdes,mat,intTab,opt);
% end
% 
%     function glb = getCouplingTerms(obj,glb,pdes,mat,intTab,opt)
%         %% interface
%         pdeDofs = obj.pdeDofs;
%         E2  = glb.feMatrices{4,1};
%         M0  = glb.feMatrices{1,3};
%         connectivity = obj.connectivity;
%         connectivityVertices = connectivity(:,1:2);                         % connectivity, only consider vertices for efficiency
%         pdeNumber = obj.pdeNumber;                                          % PDE number for each element
%         materialNumber = obj.materialNumber;                                % material number of elements
%         coordV = obj.coordV;                                                % vertex coordinates
%         elasticComponent = 2;                                               % which is the relevant dof for coupling in elastic medium
% 
%         %% check which nodes are affected
%         elasticPDEno = 0;
%         acousticPDEno = 0;
% 
%         for iPde = 1:numel(pdes)                                                    % loop PDEs
% 
%             if isa(pdes(iPde),'pdeElasticity')                                      % if PDE is elasticity
%                 elasticPDEno = iPde;                                                % remember PDE number
%             elseif isa(pdes(iPde),'pdeAcoustics')                                   % if PDE is acoustic
%                 acousticPDEno = iPde;                                               % remember PDE number
%             end
% 
%         end
% 
%         if ~elasticPDEno || ~acousticPDEno                                          % not both PDEs present
%             return                                                                  % do nothing
%         end
% 
%         dofsElastic  = pdeDofs{elasticPDEno};
%         dofsAcoustic = pdeDofs{acousticPDEno};
% 
%         nodesElastic  = ~any(isnan(dofsElastic),2);                                 % nodes belonging to elastic  PDE (logical)
%         nodesAcoustic = ~any(isnan(dofsAcoustic),2);                                % nodes belonging to acoustic PDE (logical)
%         nodesCoupled  = find(nodesElastic & nodesAcoustic);                         % nodes belonging to both
% 
%         if isa(obj,'mesh1Dcylinder')
%             isCylinder = true;
%         else
%             isCylinder = false;
%         end
% 
%         %% include coupling terms
%         for iNodes = 1:numel(nodesCoupled)                                          % loop coupled nodes
%             nodeNo = nodesCoupled(iNodes);                                          % current node number
%             eles = connectivityVertices == nodeNo;                                  % elements containing node (logical)
%             eleNos = find(any(eles,2));                                             % element numbers
%             eleElastic  = eleNos(pdeNumber(eleNos) == elasticPDEno);                % elastic element
%             eleAcoustic = eleNos(pdeNumber(eleNos) == acousticPDEno);               % acoustic element
%             yPosE = mean(coordV(connectivityVertices(eleElastic,:),2));             % midpoint of elastic element
%             yPosA = mean(coordV(connectivityVertices(eleAcoustic,:),2));            % midpoint of acoustic element
%             if yPosE > yPosA                                                        % elastic element is above acoustic element
%                 signCoupling =  1;                                                  % call this the standard sign convention
%             else                                                                    % acoustic element is above elastic element
%                 signCoupling = -1;                                                  % opposite sign
%             end
%             if isCylinder                                                   % mesh is a cylinder
%                 r = obj.coord(nodeNo,2);                                    % coupling term must be multiplied by radius at node
%             else                                                            % plate
%                 r = 1;                                                      % no extra factor
%             end
% 
%             nodeDofE = dofsElastic(nodeNo,elasticComponent);                        % relevant elastic dof of current node
%             nodeDofA = dofsAcoustic(nodeNo);                                        % relevant acoustic dof of current node
%             rho = mat(materialNumber(eleAcoustic)).parameters.rho;
%             K   = mat(materialNumber(eleAcoustic)).parameters.K;
% 
%             E2(nodeDofE,nodeDofA)= signCoupling*r;
%             M0(nodeDofA,nodeDofE)= -rho*K*signCoupling*r;
%         end
% 
%         %% output
%         glb.feMatrices{4,1} = E2;
%         glb.feMatrices{1,3} = M0;
% 
% 
%     end
end
