{ 
	if (NR==77 || NR==215) { 
			print("\\newpage"); print $0 
	} else { 
		print $0 
	}
}
