using SpidersProper

result = spiders_proper(1.25e-6, 1024; config=SpidersConfig(reference_pinhole=true))
intensity = spiders_intensity(result)

println("Final sampling: ", result.wavefront.sampling_m, " m/pixel")
println("Lyot main-pupil transmission: ", result.lyot_stop.main_pupil_transmission)
println("Lyot pinhole transmission: ", result.lyot_stop.pinhole_transmission)
println("Peak final intensity: ", maximum(intensity))
