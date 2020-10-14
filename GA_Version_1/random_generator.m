function rand_ls = random_generator(num_to_gen) % This function by default will generate random numbers between 0 - 1
    a = 0;
    b = 1;
    rand_ls = (b-a).*rand(num_to_gen,1) + a; % This is a num_to_gen * 1 column vector, we need to flip it to row vector
    rand_ls = transpose(rand_ls);
end