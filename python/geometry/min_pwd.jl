using JuMP, Pajarito, CPLEX, Mosek


function min_pwd(n, m, w, bounds, poly_bounds)
	println("do ye like dogs?")
	mip_solver_drives = true
	rel_gap = 1e-5
	mip_solver = CplexSolver(
	    CPX_PARAM_SCRIND=(mip_solver_drives ? 1 : 0),
	    # CPX_PARAM_SCRIND=1,
	    CPX_PARAM_EPINT=1e-8,
	    CPX_PARAM_EPRHS=1e-7,
	    CPX_PARAM_EPGAP=(mip_solver_drives ? 1e-5 : 1e-9)
	)

	conic_solver = MosekSolver(LOG=0)

	micp_solver = PajaritoSolver(
	    mip_solver_drives=mip_solver_drives,
	    log_level=3,
	    rel_gap=rel_gap,
		mip_solver=mip_solver,
		cont_solver=conic_solver,
	)

	mod = Model(solver = micp_solver)

	#distance variables
	@variable(mod, ds[1:binomial(n,2)-n])
	#x variables
	@variable(mod, xs[1:n,1:2])
	#eta indicator variables
	@variable(mod, etas[1:m], Bin)	

	l = 1
	for i = 1:n
		#constrain min/max values
		box = bounds[i]
		@constraint(mod, box[1:2] .<= xs[i] .<= box[3:4])
		for j = i+1:n
			#couple ds to xs and ys	
			@constraint(mod, [ds[l], xs[i,1] - xs[j,1], xs[i,2] - xs[j,2]] in SecondOrderCone())
			l = l+1
		end
	end

	#constrain points to be within region
	t = 0
	for i = 1:n
		poly_bound_i = poly_bounds[i]
		mk = len(poly_bound)
		t_start = t
		for bound in poly_bound_i
			A = bound['A']
			b = bound['b']
			@constraint(mod,A*[xs[i], etas[t]] .<= b)
			t = t+1
		end
		#constrain eta to sum to at least 1	
		@constraint(mod, sum(etas[t_start:t-1]) >= 1)
	end

	#objective
	@objective(mod, Min, w*ds)
	solve(mod)
	return getvalue(xs), getvalue(w*ds)
end


