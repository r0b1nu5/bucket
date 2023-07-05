using ReinforcementLearning, Statistics, Flux 

Base.@kwdef mutable struct CliffEnv <: AbstractEnv
	reward::Int64 = 0
	state::Vector{Int64} = [1,1]
	is_terminated::Bool = false
end

RLBase.action_space(env::CliffEnv) = 1:4
RLBase.state_space(env::CliffEnv) = 1:48
RLBase.reward(env::CliffEnv) = env.reward
RLBase.state(env::CliffEnv) = env.state
RLBase.is_terminated(env::CliffEnv) = env.is_terminated
RLBase.reset!(env::CliffEnv) = env.state,env.reward,env.is_terminated = [1,1],0,false

function (x::CliffEnv)(action)
	acts = ([1,0],[-1,0],[0,1],[0,-1])
	stts = Tuple([i,j] for i in 1:4 for j in 1:12)
	tsst = Dict{Vector{Int64},Int64}()
	for i in 1:48
		tsst[stts[i]] = i
	end

	s = stts[x.state] + acts[action]
	r = x.reward - 1
	if s[1] < 1
		s[1] = 1
	elseif s[1] > 4
		s[1] = 4
	end
	if s[2] < 1
		s[2] = 1
	elseif s[2] > 12
		s[2] = 12
	end
	if s[1] == 1 && 1 < s[2] < 12
		@info "Fell from the cliff!"
		s = [1,1]
		r -= 99
	end
	if s == [1,12]
		x.is_terminated = true
	end
	x.state = tsst[s]
	x.reward = r
end





