use rgsl::gamma_beta::gamma::gamma;

pub fn bessel_j_smallz(v: f64, z: f64) -> f64 {
    let zhalf = z / 2.0;
    let mut result = 0.0;
    let mut factorial_term = 1.0;
    let mut gamma_term = gamma(v + 1.0);
    result += 1.0 / gamma_term;
    gamma_term *= v + 1.0;
    result -= zhalf.powi(2) / gamma_term;
    for k in (2..21).step_by(2) {
        factorial_term *= k as f64;
        gamma_term *= k as f64;
        result += zhalf.powi(2 * k) / (factorial_term * gamma_term);
        factorial_term *= (k + 1) as f64;
        gamma_term *= (k + 1) as f64;
        result -= zhalf.powi(2 * k + 2) / (factorial_term * gamma_term);
    }
    result * zhalf.powf(v)
}

#[cfg(test)]
mod tests {
    use super::bessel_j_smallz;
    #[test]
    fn it_works() {
        println!("{}", bessel_j_smallz(0.0, 1.0));
        assert!((bessel_j_smallz(0.0, 1.0) - 0.765198).abs() < 1.0e-5);
    }
}
