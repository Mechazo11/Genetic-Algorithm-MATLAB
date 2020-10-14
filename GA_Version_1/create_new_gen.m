function new_gen = create_new_gen(mmx)
    % Initialize working variables
    row_idx = size(mmx,1);
    col_idx = size(mmx,2);
    cc_mat = zeros(row_idx, col_idx); % Temporary matrix to hold children in 3 by 4 shape
    ccx1 = [];
    ccx2 = [];
    % Create (N/2) random numbers
    cross_rand = random_generator(row_idx);
    
    % Select one row, fill in columns sequentially
    for i = 1:row_idx
        %fprintf("Cross over iteration %d\n", i);
        in_mmx = mmx(i,:);
        in_cross_rand = cross_rand(:,i);
        [ccx1,cxx2] = genetic_crossover(in_mmx, in_cross_rand);
        cons_mat = [ccx1, cxx2]; % Convert two row vector in 1
        % Create new_gen matrox
        for k = 1: col_idx
            cc_mat(i,k) = cons_mat(1,k);
        end
        new_gen = cc_mat; % Shape (3 by 4) for N =6 and n = 2
    end
    % Now we reshape 3 by 4 matrix into 6 by 2 to match out table
    
end

% ------------------------------------------------
function [child1, child2] = genetic_crossover(couple_row, rand_num)
    nn = 8; % We have 8 bit representation of the number
    U = 1;
    L = -1;
    J = 255;
    first_8 = 1:8;
    second_8 = 9:16;
    
    % Unpack row?
    %xxx = mmx(1,:);
    %xxx(:,3)
    
    % Unpack row vector into four gene sequences
    % mx1, mx2 --> male_parent (x1, x2)
    % fx1, fx2 --> female_parent (x1,x2)
    mx1 = couple_row(:,1);
    mx2 = couple_row(:,2);
    fx1 = couple_row(:,3);
    fx2 = couple_row(:,4);
    
    % Encode scheme
    % variable -> base_10 -> 8bit binary -> 16bit parent
    
    % Convert variable to X10
    mx1 = VartoX10(mx1,U,L,J);
    mx2 = VartoX10(mx2,U,L,J);
    fx1 = VartoX10(fx1,U,L,J);
    fx2 = VartoX10(fx2,U,L,J);
    
    % Generate a cross over point
    % 2 8bit genes makes one 16 bit geneome sequence. Thus, we can choose
    % a cross over point in (16 - 1) locations
    % CrossoverIndex=randperm(N,1)
    N= 16 - 1; 
    CrossoverIndex = cast(rand_num * N,'uint16');
    %fprintf("Cross-over point %d\n\n", CrossoverIndex);
    
    % Encode decimal to binary
    mx1 = de2bi(mx1, nn);
    mx2 = de2bi(mx2, nn);
    fx1 = de2bi(fx1, nn);
    fx2 = de2bi(fx2, nn);
    
    % Make 16 bit genome sequence (each design variable is 8bit, placed
    % side by side makes them 16 bit)
    male_gen = [mx1, mx2];
    female_gen = [fx1, fx2];
    
    % Print out 16 bit male and female gene sequence
    %disp("Male"); print_genome(male_gen)
    %disp("Female"); print_genome(female_gen)
    
    % Crossover, binary values
    child1 = [male_gen(1:CrossoverIndex) female_gen(CrossoverIndex+1:end)];
    child2 = [female_gen(1:CrossoverIndex) male_gen(CrossoverIndex+1:end)];
    
    % Print out 16 bit child gene sequence
    %disp("Child 1"); print_genome(child1);
    %disp("Child 2"); print_genome(child2);
    
    % Decode scheme
    %variable <- base_10 <- 8bit binary <- 16bit child
    [cx1, cx2, dx1, dx2] = X32toBin(child1, child2, first_8, second_8); % MATLAB's bit get only works on an integer
    
    % Check if cx1, cx2, dx1, dx2 are in binary
    %cx1;
    %cx2;
    %dx1;
    %dx2;
    
    cx1 = BintoVar(cx1, U,L,J);
    cx2 = BintoVar(cx2, U,L,J);
    dx1 = BintoVar(dx1, U,L,J);
    dx2 = BintoVar(dx2, U,L,J);
    
    % Create the 'children' row vectors
    child1 = [cx1, cx2]; 
    child2 = [dx1, dx2]; 
end

% Dependency for genetic_crossover

% ------------------------------------------------
% Xreal -- value in scale defined for the design variable
function [Xreal] = BintoVar(Xbin, U,L,J)
    Xreal = bi2de(Xbin); 
    % For some reason this has to be explicitly casted to double
    Xreal = double(Xreal);
    % Scale base-10 value to design variable scale
    Xreal = X10toVar(Xreal,U,L,J); 
end


% ------------------------------------------------
function [Xreal] = X10toVar(X10,UU,LL,JJ)
  Xreal = LL + ((UU - LL)/JJ) * X10;
  % MATLAB cannot convert a float into binary
  Xreal = double(Xreal);
end

% ------------------------------------------------
% X10 must be fininte, non-negative integer
function [X10] = VartoX10(XVar,U,L,J)
  X10 = ((XVar - L) * J)/(U - L);
  X10 = int16(X10);
end

%-----------------------------------------------
% Function to automate MATLAB's way of finding and extracting binary bits
% from an integer number
function [cx1, cx2, dx1, dx2] = X32toBin(child1, child2, first_8, second_8)
    child1 = bi2de(child1);
    child2 = bi2de(child2);
    cx1 = bitget(child1, first_8);
    cx2 = bitget(child1, second_8);
    dx1 = bitget(child2, first_8);
    dx2 = bitget(child2, second_8);
end

% ----------------------------------------------
function print_genome(P)
    fprintf('Genome : [');
    fprintf('%g ', P);
    fprintf(']\n');
    P = bi2de(P);
    %fprintf("Genome in decimal --> %d\n\n", P);
end
